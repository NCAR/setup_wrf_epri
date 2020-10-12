#!/usr/bin/env python3
'''
MIT License

Copyright (c) 2020 Jared A. Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

setup_wrf_wrapper.py

Written by: Jared A. Lee (jaredlee@ucar.edu)
Written on: 2 Oct 2020

This script is called by a crontab job with 3 arguments: WRF start hour, WRF simulation length, and IC/LBC delta hours.
This script gets the current UTC date to pair with the specified WRF start hour to provide the necessary arguments for calling setup_wrf.py.

'''

import argparse
import pathlib
import datetime as dt
import os
import sys
import time

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('init_hr', help='two-digit hour of WRF start time (e.g., 06)')
	parser.add_argument('-s', '--sim_length_h', default=6, type=int, help='integer number of hours for WRF simulation (default: 6)')
	parser.add_argument('-f', '--icbc_fc_dt', default=0, type=int, help='integer number of hours (default: 0)')

	args = parser.parse_args()
	init_hr = args.init_hr
	sim_hrs = args.sim_length_h
	icbc_fc_dt = args.icbc_fc_dt

	return init_hr, sim_hrs, icbc_fc_dt

def main(init_hr, sim_hrs, icbc_fc_dt):
	utcnow = dt.datetime.utcnow()
	date_str = utcnow.strftime('%Y%m%d')
	date_hh = date_str+'_'+init_hr

	## If this script ever moves to a different folder or system, these paths will need to be changed, as will paths in setup_wrf.py
	script_dir = pathlib.Path('/','glade','u','home','jaredlee','tools','python','setup_wrf')
	log_dir = pathlib.Path('/','glade','scratch','jaredlee','EPRI','Phase3','cron')
	log_file_run = log_dir.joinpath('cron_setup_wrf_'+date_hh+'_run.log')
	log_file_arc = log_dir.joinpath('cron_setup_wrf_'+date_hh+'_arc.log')

	## Is a loop counter needed as an additional failsafe to prevent the script from running forever?
	## Or is the status variable and checking the log file for success/error sufficient?
	wps_wrf_status = False
	while wps_wrf_status == False:
		run_cmd = str(script_dir.joinpath('setup_wrf.py'))+' '+date_hh+' -s '+str(sim_hrs)+' -f '+str(icbc_fc_dt)+' -r > '+str(log_file_run)
		arc_cmd = str(script_dir.joinpath('setup_wrf.py'))+' '+date_hh+' -s '+str(sim_hrs)+' -f '+str(icbc_fc_dt)+' -a > '+str(log_file_arc)
		print('Executing command:\n'+run_cmd)
		os.system(run_cmd)

		## Are the icbc files not being found on Google Cloud or RDA? If so, then try again with icbc_fc_dt + 1
		## Note that for HRRR on Google Cloud, while most files start arriving ~50 min after initialization, some don't arrive for 1h 25m on occasion
		if 'ERROR: HTTP Error 404: Not Found' in open(str(log_file_run)).read():
			print('One or more of the IC/BC files not found. Trying an earlier cycle.')
			icbc_fc_dt = icbc_fc_dt + 1
			print('Setting icbc_fc_dt = '+str(icbc_fc_dt)+' h.')

		## Some other error
		elif 'ERROR' in open(str(log_file_run)).read():
			print('There was an error running setup_wrf.py. Check '+str(log_file_run)+' for details.')
			print('Archive step not run. Exiting script.')
			sys.exit()

		## It worked!
		elif 'SUCCESS! WRF completed successfully.' in open(str(log_file_run)).read():
			wps_wrf_status = True

		## Something else happened. Check log file, and add new elif branch here if needed.
		else:
			print('Something unknown happened while running setup_wrf.py. Check '+str(log_file_run)+' for details. Exiting script.')
			sys.exit()

	## Assuming everything ran correctly, now run the archiving step, after giving the file system opportunity to collect itself
	time.sleep(1)
	print('Executing command:\n'+arc_cmd)
	os.system(arc_cmd)

	print('setup_wrf_master.py has completed.')


if __name__ == '__main__':
	init_hr, sim_hrs, icbc_fc_dt = parse_args()
	main(init_hr, sim_hrs, icbc_fc_dt)
