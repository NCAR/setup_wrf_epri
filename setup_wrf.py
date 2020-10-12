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

setup_wrf.py

Written by: Jared A. Lee
            jaredlee@ucar.edu
Written on: 20 Sep 2017

This script links to (and downloads, if necessary) all the files needed to run WPS/WRF.
Each program in the WPS/WRF workflow can be optionally executed if its option is set to True.
WRF output files can also be optionally moved to an archival directory (arc_dir).
This script is built on my original setup_wrf.bash, which is clunky for text replacement with sed/awk.

Edits:
======
27 Jun 2018
-- Provided option for ungribbing two datasets, one for general IC/BC and one for soil fields (EPRI/Phase2).
27 Jul 2018
-- Created command-line argument handling for variables 'beg_dt' and 'sim_hrs'.
   Variable end_dt now calculated automatically from 'beg_dt' and 'sim_hrs'.
	Added option for compressing archived files with gzip automatically.
29 May 2019
-- Replaced numerous os.system, os.mkdir, os.remove, os.path.* calls with:
	- pathlib.Path objects/functions
	- shutil.move/shutil.copy
	Some os.system calls remain to execute external scripts. May replace with subprocess eventually.
08 Nov 2019
-- Added ability to copy GEOS-5 forecast files from a remote server, and run a Fortran executable to prepare
   auxinput files for WRF, to allow GEOS-5 forecast AOD to be imposed on the WRF simulation.
03 Mar 2020
-- Added looping and options to handle multiple initialization times (with integer hour stride).
   Added some additional command-line optional arguments.
	Put parse_args into its own function, the core stuff into a main function, & added __main__ block at end.
	Added run time reporting in __main__ block at the end.
01 Apr 2020
-- Added necessary options for running WRF for EPRI/Phase3 project.
02 Oct 2020
-- Added realtime flag to keep the script running until WRF is successful or terminates with an error.
-- Added icbc_fc_dt as an optional command line argument to facilitate rerunning the script with an older IC/LBC cycle.
08 Oct 2020
-- Cleaned up old date format lines, replaced string position extraction with pd.to_datetime.
-- Cleaned up some obsolete commented code chunks.
'''

import os
import pathlib
import shutil
import datetime as dt
import pytz
import numpy as np
import pandas as pd
import sys
import fileinput
import time
import argparse
import wget

def parse_args():
	## Parse the command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('init_dt_first', help='beginning date/time of first WRF simulation [YYYYMMDD_HH]')
	parser.add_argument('-l', '--init_dt_last', default=None, help='beginning date/time of last WRF simulation [YYYYMMDD_HH] (default: same as init_dt_first)')
	parser.add_argument('-i', '--init_stride_h', default=24, type=int, help='integer number of hours between forecast cycles (default: 24)')
	parser.add_argument('-s', '--sim_length_h', default=48, type=int, help='integer number of hours for WRF simulation (default: 48)')
	parser.add_argument('-f', '--icbc_fc_dt', default=0, type=int, help='integer number of hours prior to WRF init time for IC/LBC model cycle (default: 0)')
	parser.add_argument('-r', '--realtime', help='flag when running in real-time to keep this script running until WRF is done', action='store_true')
	parser.add_argument('-a', '--archive', help='flag to archive wrfout, wrfinput, wrfbdy, and namelist files off of scratch', action='store_true')

	args = parser.parse_args()
	init_dt_first = args.init_dt_first
	init_dt_last  = args.init_dt_last
	init_stride_h = args.init_stride_h
	sim_hrs = args.sim_length_h
	icbc_fc_dt = args.icbc_fc_dt

	if len(init_dt_first) != 11:
		print('ERROR! Incorrect format for positional argument init_dt_first. Exiting!')
		parser.print_help()
		sys.exit()
	elif init_dt_first[8] != '_':
		print('ERROR! Incorrect format for positional argument init_dt_first. Exiting!')
		parser.print_help()
		sys.exit()

	if init_dt_last != None:
		if len(init_dt_last) != 11:
			print('ERROR! Incorrect length for positional argument init_dt_last. Exiting!')
			parser.print_help()
			sys.exit()
		elif init_dt_last[8] != '_':
			print('ERROR! Incorrect format for positional argument init_dt_last. Exiting!')
			parser.print_help()
			sys.exit()
	else:
		init_dt_last = init_dt_first

	rt_flag = False
	if args.realtime:
		rt_flag = True

	arc_flag = False
	if args.archive:
		arc_flag = True

	return init_dt_first, init_dt_last, init_stride_h, sim_hrs, icbc_fc_dt, rt_flag, arc_flag

def wget_error(error_msg, now_time_beg):
	print('ERROR: '+error_msg)
	print('Check if an earlier cycle has the required files and adjust icbc_fc_dt if necessary. Exiting!')
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('')
	print('Script completed.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot))
	print('')
	sys.exit()

def main(init_dt_first, init_dt_last, init_stride_h, sim_hrs, icbc_fcst_dt, rt_flag, arc_flag, now_time_beg):
	project = 'EPRI/Phase3'

	del_run_dir		= False	# Delete existing run directory?
	link_wps			= True	# Link to WPS executables?
	link_wrf			= True	# Link to WRF files?
	link_met_em		= False	# Link to pre-existing metgrid output data (met_em*)
	link_realwrf	= False	# Link to real & wrf from experiment-specific directory?
	link_joinwrf	= False	# Link to joinwrf?
	link_wrfinbdy	= False	# Link to pre-existing wrfinput and wrfbdy files?
	link_wrflowinp	= False	# Link to pre-existing wrflowinp files?
	write_ts			= False	# Write time series output files and link to tslist file for time series locations?
	is_fdda			= False	# Community FDDA run?
	is_rtfdda		= False	# RT-FDDA run?
	is_skebs_ens	= False	# SKEBS ensemble run?
	get_icbc			= False	# Download/link to IC/BC grib data?
	icbc_anal		= False	# Use analysis files for IC/BC?
	icbc_fcst		= False	# Use forecast files for IC/BC?
	run_geogrid		= False	# Run geogrid for this case?
	run_ungrib		= False	# Run ungrib for this case?
	ungrib_soil		= False	# Run ungrib separately for soil data?
	run_avg_tsfc	= False	# Run avg_tsfc for this case?
	run_metgrid		= False	# Run metgrid for this case?
	sub_real			= False	# Submit real for this case?
	aero_scp_data	= False	# scp aerosol data from a remote server?
	impose_aerosol	= False	# Impose aerosols?
	sub_wrf			= False	# Submit wrf for this case?
	arc_wrf			= False	# Archive wrfout, wrfinput, wrfbdy files for this case?
	arc_wrfout		= True	# Archive wrfout files? May not want to if processed through UPP
	arc_uppout		= False	# Archive UPP output grib2 files?
	arc_gzip			= False	# Compress archived files with gzip

	lead_stride = 1	# default stride in lead hours for downloading IC/BC files if available (may only want every 3 h from NOMADS for GFS)

	time_fmt_hpss_date_dir  = '%Y%m%d'
	time_fmt_exp_dir        = '%Y-%m-%d_%H'
	time_fmt_yyyymmdd       = '%Y%m%d'
	time_fmt_yyyymmddhh     = '%Y%m%d%H'
	time_fmt_yyyymmdd_hh    = '%Y%m%d_%H'
	time_fmt_yyyymmdd_hhmm  = '%Y%m%d_%H%M'

	## Date/time manipulation
	init_datetime_first = pd.to_datetime(init_dt_first, format=time_fmt_yyyymmdd_hh)
	init_datetime_last  = pd.to_datetime(init_dt_last,  format=time_fmt_yyyymmdd_hh)

	dif_datetime = init_datetime_last - init_datetime_first
	dif_datetime_h = dif_datetime.days * 24 + dif_datetime.seconds // 3600

	n_init_times = dif_datetime_h // init_stride_h + 1  # / = floating point division, // = integer (floor) division

	this_datetime = init_datetime_first
	init_datetime_all = ['' for x in range(n_init_times)]
	init_yrs = ['' for x in range(n_init_times)]
	init_mos = ['' for x in range(n_init_times)]
	init_dys = ['' for x in range(n_init_times)]
	init_hrs = ['' for x in range(n_init_times)]

	for tt in range(n_init_times):
		init_datetime_all[tt] = this_datetime.replace(tzinfo=pytz.utc)
		init_yrs[tt] = this_datetime.strftime('%Y')
		init_mos[tt] = this_datetime.strftime('%m')
		init_dys[tt] = this_datetime.strftime('%d')
		init_hrs[tt] = this_datetime.strftime('%H')

		this_datetime = this_datetime + dt.timedelta(hours=init_stride_h)

	## Loop over forecast cycles/initializations
	for ii in range(n_init_times):
		init_datetime = init_datetime_all[ii]
		init_yr = init_yrs[ii]
		init_mo = init_mos[ii]
		init_dy = init_dys[ii]
		init_hr = init_hrs[ii]
		beg_dt = init_yr+'-'+init_mo+'-'+init_dy+'_'+init_hr

		sub_dir		= pathlib.Path('/','glade','u','home','jaredlee','tools','submit_scripts') # LSF/PBS submit scripts
		joinwrf_dir	= pathlib.Path('/','glade','u','home','jaredlee','tools','joinwrf')  # Dir containing joinwrf exe and perl script (for tiled WRF runs)

		if project == 'EPRI/Phase3':
			num_dom		= 1
			ens_mem		= 1	# number of ensemble members
			if sim_hrs <= 6:
				is_nowcast = True
				is_dayahead = False
			elif sim_hrs == 42:
				is_nowcast = False
				is_dayahead = True
			else:
				print('ERROR: Please provide "-s 6" or "-s 42" as an optional argument to setup_wrf.py.')
				print('       -s 6 = nowcast run; -s 42 = dayahead run')
				print('Also note that changing this from 42 will require changing default WPS and WRF namelists.')
				print('Exiting!')
				sys.exit()
			if is_nowcast and not is_dayahead:
				run_dir	= pathlib.Path('/','glade','scratch','jaredlee',project,'wrf','nowcast',beg_dt)
				arc_dir	= pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','nowcast',beg_dt)
				nml_dir	= pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','nowcast','config','nml_templates')
				sub_dir  = pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','nowcast','config') 
				realwrf_dir	= pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','nowcast','config')
				icbc_model	= 'HRRR'
			elif is_dayahead and not is_nowcast:
				run_dir	= pathlib.Path('/','glade','scratch','jaredlee',project,'wrf','dayahead',beg_dt)
				arc_dir	= pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','dayahead',beg_dt)
				nml_dir	= pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','dayahead','config','nml_templates')
				sub_dir	= pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','dayahead','config')
				realwrf_dir = pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','dayahead','config')
				icbc_model	= 'GFS'
#				gfs_source  = 'RDA'
				gfs_source  = 'NOMADS'	# only use this for near-real-time runs, as NOMADS only retains last 10 days
			else:
				print('ERROR: is_nowcast and is_dayahead have the same value. Exiting!')
				sys.exit()
			tslist_file	= realwrf_dir.joinpath('tslist')
			wrf_dir		= pathlib.Path('/','glade','u','home','jaredlee','programs','WRF_v4.2','run')
			wps_dir		= pathlib.Path('/','glade','u','home','jaredlee','programs','WPS_v4.2')
			geo_dir		= pathlib.Path('/','glade','p','ral','wsap','jaredlee',project,'wrf','config')
			if icbc_model == 'HRRR':
#				icbc_fc_dt	= 2	# time delta (hours) for pulling IC/BC files from a cycle prior to beg_dt, if desired
				icbc_dir		= pathlib.Path('https://storage.googleapis.com', 'high-resolution-rapid-refresh')	# Google public cloud storage (free)
				soil_dir		= icbc_dir
				icbc_vtable	= 'Vtable.RAP-HRRR_hybrid'
				soil_vtable	= 'Vtable.RAP-HRRR_soil_only'
				icbc_nml 	= 'namelist.wps.hrrr.hybr'
				soil_nml 	= 'namelist.wps.hrrr.soil'
				ungrib_soil	= True
			elif icbc_model == 'GFS':
#				icbc_fc_dt	= 0
				if gfs_source == 'RDA':
					## NOTE: GLADE RDA does not get a given day's GFS until 6am-7am MDT (12-13 UTC) the next day
					icbc_dir 	= pathlib.Path('/','glade','collections','rda','data','ds084.1')	# glade RDA (0.25-deg, 3-hourly GFS)
				elif gfs_source == 'NOMADS':
					## NOTE: NCEP NOMADS only retains 9 days of model runs + current day. 42 h of 06z GFS in by 0938z.
					icbc_dir		= pathlib.Path('https://nomads.ncep.noaa.gov','pub','data','nccf','com','gfs','prod')	# NCEP NOMADS
					lead_stride = 3	# change this to 1 if it's desired to download every hour of GFS from NOMADS for IC/BC files
				icbc_vtable	= 'Vtable.GFS'
				icbc_nml 	= 'namelist.wps.gfs'
				ungrib_soil	= False
			grib_dir		= pathlib.Path('/','glade','scratch','jaredlee',project,'wps','grib',icbc_model.lower(),beg_dt)
			ungrib_dir	= pathlib.Path('/','glade','scratch','jaredlee',project,'wps','ungrib',beg_dt)
			met_dir		= pathlib.Path('/','glade','scratch','jaredlee',project,'wps','metgrid',beg_dt)

			write_ts		= False
			get_icbc		= True
			icbc_fcst	= True
			run_geogrid	= False
			run_ungrib	= True
			run_metgrid	= True
			sub_real		= True
			aero_scp_data	= False	# will need to look at this for CAMS vs GEOS-5
			impose_aerosol	= False
			sub_wrf		= True
		else:
			print('ERROR: No options set for project = '+project+'. Exiting!')
			sys.exit()

		if arc_wrf == True:
			arc_flag	= True

		## If the archive flag was included in the command line, then just do that
		if arc_flag:
			link_wps			= False
			link_wrf			= False
			link_met_em		= False
			link_realwrf	= False
			link_joinwrf	= False
			link_wrfinbdy	= False
			link_wrflowinp	= False
			get_icbc			= False
			icbc_anal		= False
			icbc_fcst		= False
			run_geogrid		= False
			run_ungrib		= False
			run_avg_tsfc	= False
			run_metgrid		= False
			sub_real			= False
			impose_aerosol	= False
			sub_wrf			= False
			arc_wrf			= True

		## Parse the date strings into their components
		beg_yr = init_yr
		beg_mo = init_mo
		beg_dy = init_dy
		beg_hr = init_hr

		## Date/time manipulation
		beg_datetime	= init_datetime
		end_datetime	= beg_datetime + dt.timedelta(hours=sim_hrs)
		beg_date = beg_datetime.date()
		beg_time = beg_datetime.time()
		end_date	= end_datetime.date()
		end_time	= end_datetime.time()
		end_dt	= end_datetime.strftime(time_fmt_exp_dir)

		beg_dt_yyyymmddhh	= beg_datetime.strftime(time_fmt_yyyymmddhh)
		end_dt_yyyymmddhh = end_datetime.strftime(time_fmt_yyyymmddhh)
		beg_dt_yyyymmdd_hhmm = beg_datetime.strftime(time_fmt_yyyymmdd_hhmm)
		end_dt_yyyymmdd_hhmm = end_datetime.strftime(time_fmt_yyyymmdd_hhmm)

		end_yr	= end_dt[0:4]
		end_mo	= end_dt[5:7]
		end_dy	= end_dt[8:10]
		end_hr	= end_dt[11:]

		if get_icbc and icbc_fcst:
			icbc_fc_datetime	= beg_datetime - dt.timedelta(hours=icbc_fc_dt)
			icbc_fc_cyc	= icbc_fc_datetime.strftime(time_fmt_exp_dir)

		exp_duration	= end_datetime - beg_datetime
		exp_duration_h	= (exp_duration.seconds // 3600) + (exp_duration.days * 24)	# this is in UTC, so don't need to worry about DST changes

		## If there's a separate end date for ungrib/avg_tsfc, then do some date/time manipulation on that
		try:
			end_ungrib
		except NameError:
			end_ungrib_exists = False
		else:
			end_ungrib_exists = True
			end_yr_u = end_ungrib[0:4]
			end_mo_u = end_ungrib[5:7]
			end_dy_u = end_ungrib[8:10]
			end_hr_u = end_ungrib[11:]
			end_date_u		= dt.date(int(end_yr_u), int(end_mo_u), int(end_dy_u))
			end_datetime_u	= dt.datetime(int(end_yr_u), int(end_mo_u), int(end_dy_u), int(end_hr_u), 0, 0)
			if end_datetime_u < end_datetime:
				if (run_ungrib or run_avg_tsfc):
					print('ERROR: Chosen end date/time for ungrib data ('+end_ungrib+') precedes the simulation end date ('+end_dt+').')
					print('       This will result in insufficient IC/BC model data being processed by ungrib or avg_tsfc. Exiting!')
					sys.exit()
				else:
					print('WARNING: Chosen end date/time for ungrib data ('+end_ungrib+') precedes the simulation end date ('+end_dt+').')
					print('         Change this setting if you want to run ungrib and/or avg_tsfc.')

		## Does the run_dir exist? If not, then create it.
		## But first, do we want to blitz the run directory and start over?
		if del_run_dir:
			if pathlib.Path.is_dir(run_dir):
				print('Removing old run directory '+str(run_dir))
			shutil.rmtree(run_dir)
		pathlib.Path(run_dir).mkdir(parents=True, exist_ok=True)

		## Does the arc_dir and its subdirectories exist? If not, then create it.
		pathlib.Path(arc_dir).mkdir(parents=True, exist_ok=True)
		arc_dir.joinpath('config').mkdir(parents=True, exist_ok=True)
		arc_dir.joinpath('wrfout').mkdir(parents=True, exist_ok=True)
		arc_dir.joinpath('plots').mkdir(parents=True, exist_ok=True)
		arc_dir.joinpath('uppout').mkdir(parents=True, exist_ok=True)

		## Go to run_dir
		os.chdir(str(run_dir))

		## Link to WPS executables
		if link_wps:
			print('Linking WPS files in '+str(run_dir)+' to '+str(wps_dir))
			if pathlib.Path('geogrid').is_symlink():
				pathlib.Path('geogrid').unlink()
			if pathlib.Path('ungrib').is_symlink():
				pathlib.Path('ungrib').unlink()
			if pathlib.Path('metgrid').is_symlink():
				pathlib.Path('metgrid').unlink()
			if pathlib.Path('util').is_symlink():
				pathlib.Path('util').unlink()
			if pathlib.Path('geogrid.exe').is_symlink():
				pathlib.Path('geogrid.exe').unlink()
			if pathlib.Path('ungrib.exe').is_symlink():
				pathlib.Path('ungrib.exe').unlink()
			if pathlib.Path('metgrid.exe').is_symlink():
				pathlib.Path('metgrid.exe').unlink()
			if pathlib.Path('avg_tsfc.exe').is_symlink():
				pathlib.Path('avg_tsfc.exe').unlink()
			if pathlib.Path('link_grib.csh').is_symlink():
				pathlib.Path('link_grib.csh').unlink()

			for x in pathlib.Path(wps_dir).glob('*.exe'):
				pathlib.Path(x.name).symlink_to(x)
			pathlib.Path('geogrid').symlink_to(wps_dir.joinpath('geogrid'))
			pathlib.Path('ungrib').symlink_to(wps_dir.joinpath('ungrib'))
			pathlib.Path('metgrid').symlink_to(wps_dir.joinpath('metgrid'))
			pathlib.Path('util').symlink_to(wps_dir.joinpath('util'))
			pathlib.Path('avg_tsfc.exe').symlink_to(wps_dir.joinpath('util','avg_tsfc.exe'))
			pathlib.Path('link_grib.csh').symlink_to(wps_dir.joinpath('link_grib.csh'))

		## Link to WRF files. Link to everything in wrf_dir/run, but delete link to namelist to prepare for new namelist being copied over
		if link_wrf:
			print('Linking WRF files in '+str(run_dir)+' to '+str(wrf_dir))
#			for fname in glob.glob(wrf_dir+'/*'):	# Returns the file names matching the pattern with the absolute path prepended
#			for fname in os.listdir(wrf_dir):		# Returns just the file names within the directory
#				symlink_force(wrf_dir+'/'+fname, fname)
			for x in pathlib.Path(wrf_dir).glob('*'):
				if not pathlib.Path(x.name).is_symlink():
					pathlib.Path(x.name).symlink_to(x)
			if pathlib.Path('namelist.input').is_symlink():
				pathlib.Path('namelist.input').unlink()
			run_dir.joinpath('rsl_real').mkdir(parents=True, exist_ok=True)
			run_dir.joinpath('rsl_wrf').mkdir(parents=True, exist_ok=True)
			try: vars_io
			except NameError: vars_io = None
			if vars_io != None:
				if vars_io.is_file():
					shutil.copy(str(vars_io), 'vars_io.txt')


		## Link to experiment-specific real & wrf executables?
		if link_realwrf:
			print('Linking experiment-specific real.exe and wrf.exe executables to '+str(run_dir))
			pathlib.Path('real.exe').unlink()
			pathlib.Path('wrf.exe').unlink()
			pathlib.Path('real.exe').symlink_to(realwrf_dir.joinpath('real.exe'))
			pathlib.Path('wrf.exe').symlink_to(realwrf_dir.joinpath('wrf.exe'))

		## Link to metgrid output
		if link_met_em:
			print('Linking metgrid output to '+str(run_dir))
			for x in pathlib.Path(metdir).glob('met_em*'):
				if not pathlib.Path(x.name).is_symlink():
					pathlib.Path(x.name).symlink_to(x)

		## Link to time series site list file
		if write_ts and not arc_flag:
			print('Linking tslist file to '+str(run_dir))
			if not pathlib.Path(tslist_file.name).is_symlink():
				pathlib.Path(tslist_file.name).symlink_to(tslist_file)
			arc_dir.joinpath('ts').mkdir(parents=True, exist_ok=True)

		## Copy PBS submit scripts
		shutil.copy(str(sub_dir.joinpath('pbs_submit_wrf.bash')), '.')
		shutil.copy(str(sub_dir.joinpath('pbs_submit_real.bash')), '.')

		## Create new start and end date namelist strings according to number of domains
		n = 1
		wps_nml_start_date_line	= " start_date = '"+beg_dt+":00:00',"
		wps_nml_end_date_line	= " end_date   = '"+end_dt+":00:00',"

		## Set default prefix and fg_name lines
		wps_nml_prefix_line  = " prefix = '"+str(ungrib_dir)+'/'+icbc_model+"',"
		wps_nml_fg_name_line = " fg_name = '"+str(ungrib_dir)+'/'+icbc_model+"',"

		## Need special handling for initializing from HRRR hybrid level data (as soil data is in pressure level files, need to run ungrib twice)
		if icbc_model == 'HRRR':
			wps_nml_prefix_line			= " prefix = '"+str(ungrib_dir)+'/'+icbc_model+"_hybr',"
			wps_nml_prefix_line_soil	= " prefix = '"+str(ungrib_dir)+'/'+icbc_model+"_soil',"
			wps_nml_fg_name_line			= " fg_name = '"+str(ungrib_dir)+'/'+icbc_model+"_hybr',"
			wps_nml_fg_name_line_soil	= " fg_name = '"+str(ungrib_dir)+'/'+icbc_model+"_hybr','"+str(ungrib_dir)+'/'+icbc_model+"_soil',"

		wps_nml_metgrid_out_line	= " opt_output_from_metgrid_path = '"+str(met_dir)+"',"
		wrf_nml_start_year_line		= ' start_year                          = '+beg_yr+','
		wrf_nml_start_month_line	= ' start_month                         = '+beg_mo+','
		wrf_nml_start_day_line		= ' start_day                           = '+beg_dy+','
		wrf_nml_start_hour_line		= ' start_hour                          = '+beg_hr+','
		wrf_nml_end_year_line		= ' end_year                            = '+end_yr+','
		wrf_nml_end_month_line		= ' end_month                           = '+end_mo+','
		wrf_nml_end_day_line			= ' end_day                             = '+end_dy+','
		wrf_nml_end_hour_line      = ' end_hour                            = '+end_hr+','
		wrf_nml_run_hours_line		= ' run_hours                           = '+str(exp_duration_h)+','
		while n < num_dom:
			wps_nml_start_date_line = wps_nml_start_date_line +" '"+beg_dt+":00:00',"
			wps_nml_end_date_line	= wps_nml_end_date_line   +" '"+beg_dt+":00:00',"
			wrf_nml_start_year_line	= wrf_nml_start_year_line +'  '  +beg_yr+','
			wrf_nml_start_month_line= wrf_nml_start_month_line+'    '+beg_mo+','
			wrf_nml_start_day_line	= wrf_nml_start_day_line  +'    '+beg_dy+','
			wrf_nml_start_hour_line	= wrf_nml_start_hour_line +'    '+beg_hr+','
			wrf_nml_end_year_line	= wrf_nml_end_year_line   +'  '  +end_yr+','
			wrf_nml_end_month_line	= wrf_nml_end_month_line  +'    '+end_mo+','
			wrf_nml_end_day_line		= wrf_nml_end_day_line    +'    '+end_dy+','
			wrf_nml_end_hour_line	= wrf_nml_end_hour_line   +'    '+end_hr+','
			n = n + 1

		## WPS namelist handling
		if link_wps:
			print('Editing dates and linking WPS namelist to '+str(arc_dir)+'/config/'+icbc_nml)
			with open(str(nml_dir.joinpath(icbc_nml)), 'r') as input_file, open(str(nml_dir.joinpath(icbc_nml+'.new')), 'w') as output_file:
				for line in input_file:
					if line.strip()[0:10] == 'start_date':
						output_file.write(wps_nml_start_date_line+'\n')
					elif line.strip()[0:8] == 'end_date':
						output_file.write(wps_nml_end_date_line+'\n')
					elif line.strip()[0:6] == 'prefix':
						output_file.write(wps_nml_prefix_line+'\n')
					elif line.strip()[0:7] == 'fg_name':
						if ungrib_soil:
							output_file.write(wps_nml_fg_name_line_soil+'\n')
						else:
							output_file.write(wps_nml_fg_name_line+'\n')
					elif line.strip()[0:28] == 'opt_output_from_metgrid_path':
						output_file.write(wps_nml_metgrid_out_line+'\n')
					else:
						output_file.write(line)

			shutil.move(str(nml_dir.joinpath(icbc_nml+'.new')), str(arc_dir.joinpath('config',icbc_nml)))
			## Always force linking to the newly updated namelist in case something changed or went wrong previously
			if pathlib.Path(icbc_nml).is_symlink():
				pathlib.Path(icbc_nml).unlink()
			pathlib.Path(icbc_nml).symlink_to(arc_dir.joinpath('config',icbc_nml))

			if ungrib_soil:
				print('Editing dates and linking WPS namelist to '+str(arc_dir)+'/config/'+soil_nml)
				with open(str(nml_dir.joinpath(soil_nml)), 'r') as input_file, open(str(nml_dir.joinpath(soil_nml+'.new')), 'w') as output_file:
					for line in input_file:
						if line.strip()[0:10] == 'start_date':
							output_file.write(wps_nml_start_date_line+'\n')
						elif line.strip()[0:8] == 'end_date':
							output_file.write(wps_nml_end_date_line+'\n')
						elif line.strip()[0:6] == 'prefix':
							output_file.write(wps_nml_prefix_line_soil+'\n')
						elif line.strip()[0:7] == 'fg_name':
							output_file.write(wps_nml_fg_name_line_soil+'\n')
						elif line.strip()[0:28] == 'opt_output_from_metgrid_path':
							output_file.write(wps_nml_metgrid_out_line+'\n')
						else:
							output_file.write(line)

					shutil.move(str(nml_dir.joinpath(soil_nml+'.new')), str(arc_dir.joinpath('config',soil_nml)))
					## Always force linking to the newly updated namelist in case something changed or went wrong previously
					if pathlib.Path(soil_nml).is_symlink():
						pathlib.Path(soil_nml).unlink()
					pathlib.Path(soil_nml).symlink_to(arc_dir.joinpath('config',soil_nml))

		## WRF namelist handling
		if link_wrf:
			print('Editing dates and linking WRF namelist to '+str(arc_dir)+'/config/namelist.input')
			if is_fdda:
				print('-- Using FDDA WRF namelist')
				wrf_nml_template = nml_dir.joinpath('namelist.input.fdda')
			elif is_rtfdda:
				print('-- Using RT-FDDA WRF namelist')
				wrf_nml_template = nml_dir.joinpath('namelist.input.rtfdda')
			else:
				wrf_nml_template = nml_dir.joinpath('namelist.input')

			with open(str(wrf_nml_template), 'r') as input_file, open(str(nml_dir.joinpath('namelist.input.new')), 'w') as output_file:
				for line in input_file:
					if line.strip()[0:9] == 'run_hours':
						output_file.write(wrf_nml_run_hours_line+'\n')
					elif line.strip()[0:10] == 'start_year':
						output_file.write(wrf_nml_start_year_line+'\n')
					elif line.strip()[0:11] == 'start_month':
						output_file.write(wrf_nml_start_month_line+'\n')
					elif line.strip()[0:9]  == 'start_day':
						output_file.write(wrf_nml_start_day_line+'\n')
					elif line.strip()[0:10] == 'start_hour':
						output_file.write(wrf_nml_start_hour_line+'\n')
					elif line.strip()[0:8]  == 'end_year':
						output_file.write(wrf_nml_end_year_line+'\n')
					elif line.strip()[0:9]  == 'end_month':
						output_file.write(wrf_nml_end_month_line+'\n')
					elif line.strip()[0:7]  == 'end_day':
						output_file.write(wrf_nml_end_day_line+'\n')
					elif line.strip()[0:8]  == 'end_hour':
						output_file.write(wrf_nml_end_hour_line+'\n')
					else:
						output_file.write(line)

			shutil.move(str(nml_dir.joinpath('namelist.input.new')), str(arc_dir.joinpath('config','namelist.input')))
			## Always force linking to the newly updated namelist in case something changed or went wrong previously
			if pathlib.Path('namelist.input').is_symlink():
				pathlib.Path('namelist.input').unlink()
			pathlib.Path('namelist.input').symlink_to(arc_dir.joinpath('config','namelist.input'))

		## If this is an ensemble (SKEBS) experiment, then set up the ensemble member directories
		if is_skebs_ens:
			## First, check that ens_mem has been set
			try:
				ens_mem
			except NameError:
				ens_mem_exists = False
				print('ERROR: is_skebs_ens = True but ens_mem is not set. Exiting!')
				sys.exit()
			else:
				ens_mem_exists = True	# If this line executes, then everything's fine
				if ens_mem >= 10000:
					print('ERROR: ens_mem >= 10000. This exceeds max ensemble size for this script currently. Exiting!')
					sys.exit()

			## Create numeric and string arrays of ensemble member numbers (with appropriate leading zeroes)
			ens_mem_arr			= np.arange(1, ens_mem+1)
			ens_mem_arr_str	= []
			for mem in ens_mem_arr:
				if ens_mem < 100:
					ens_mem_arr_str.append(str(mem).zfill(2))
				elif ens_mem < 1000:
					ens_mem_arr_str.append(str(mem).zfill(3))
				else:
					ens_mem_arr_str.append(str(mem).zfill(4))

			## Create and fill ensemble member directories, both in arc_dir and run_dir
			m = 0
			for member in ens_mem_arr_str:
				run_dir.joinpath('mem'+member).mkdir(parents=True, exist_ok=True)
				run_dir.joinpath('mem'+member,'rsl_wrf').mkdir(parents=True, exist_ok=True)
				arc_dir.joinpath('mem'+member).mkdir(parents=True, exist_ok=True)
				arc_dir.joinpath('mem'+member,'config').mkdir(parents=True, exist_ok=True)
				arc_dir.joinpath('mem'+member,'wrfout').mkdir(parents=True, exist_ok=True)
				arc_dir.joinpath('mem'+member,'plots').mkdir(parents=True, exist_ok=True)

				if link_wrf:
					print('Creating and populating run directory for SKEBS ensemble member '+member+' in '+str(run_dir)+'/mem'+member)
					os.chdir(str(run_dir)+'/mem'+member)
					for x in pathlib.Path(wrf_dir).glob('*'):
						if not pathlib.Path(x.name).is_symlink():
							pathlib.Path(x.name).symlink_to(x)
					pathlib.Path('namelist.input').unlink()
					shutil.copy('../pbs_submit_wrf.bash', '.')

					## Edit the WRF namelist to use the right SKEBS seed for this member
					with open('../namelist.input','r') as input_file, open(str(arc_dir.joinpath('mem'+member,'config','namelist.input')),'w') as output_file:
						for line in input_file:
							if line.strip()[0:4] == 'nens':
								output_file.write(' nens                                = '+str(ens_mem_arr[m])+',\n')
							else:
								output_file.write(line)

					if not pathlib.Path('namelist.input').is_symlink():
						pathlib.Path('namelist.input').symlink_to(arc_dir.joinpath('mem'+member,'config','namelist.input'))
					os.chdir(str(run_dir))
				m = m+1

		## Copy joinwrf
		if link_joinwrf:
			print('Linking joinwrf executable and copying perl script to '+str(run_dir))
			if not pathlib.Path('joinwrf').is_symlink():
				pathlib.Path('joinwrf').symlink_to(joinwrf_dir.joinpath('joinwrf'))
			shutil.copy(str(joinwrf_dir.joinpath('join_wrf.pl')), '.')
			run_dir.joinpath('indiv').mkdir(parents=True, exist_ok=True)
			run_dir.joinpath('rsl_joinwrf').mkdir(parents=True, exist_ok=True)
			run_dir.joinpath('logs_joinwrf').mkdir(parents=True, exist_ok=True)

		## Link to wrfinput and wrfbdy files?
		if link_wrfinbdy:
			print('Linking wrfinput and wrfbdy files to '+str(run_dir))
			for x in pathlib.Path(wrfinbdy_dir).glob('wrfinput*'):
				if pathlib.Path(x).is_file():
					pathlib.Path(x).unlink()
				pathlib.Path(x.name).symlink_to(x)
			if pathlib.Path('wrfbdy_d01').is_file():
				pathlib.Path('wrfbdy_d01').unlink()
			pathlib.Path('wrfbdy_d01').symlink_to(wrfinbdy_dir.joinpath('wrfbdy_d01'))

		## Link to wrflowinp files?
		if link_wrflowinp:
			print('Linking wrflowinp files to '+str(run_dir))
			for x in pathlib.Path(wrfinbdy_dir).glob('wrflowinp*'):
				if pathlib.Path(x).is_file():
					pathlib.Path(x).unlink()
				pathlib.Path(x.name).symlink_to(x)

		## Run geogrid
		if run_geogrid:
			print('Running geogrid . . .')
			if not pathlib.Path('namelist.wps').is_file():
				if pathlib.Path('namelist.wps.hybr').is_file():
					pathlib.Path('namelist.wps').symlink_to('namelist.wps.hybr')
	
			os.system('./geogrid.exe >& geogrid.out')
			if 'Successful completion' in open('geogrid.log').read():
				print('   SUCCESS! geogrid completed successfully.')
			else:
				print('   ERROR: geogrid.exe failed.')
				print('   Consult '+str(run_dir)+'/geogrid.log and geogrid.out for potential error messages.')
				print('   Exiting!')
				sys.exit()

		## Download IC/BC grib files from GLADE RDA or external website
		if get_icbc:
			grib_dir.mkdir(parents=True, exist_ok=True)
			print('Downloading IC/BC grib files to '+str(grib_dir))
			os.chdir(str(grib_dir))
			if icbc_anal and icbc_fcst:
				print('   ERROR: Both icbc_anal and icbc_fcst = True. Set one to False and re-run this script. Exiting!')
				sys.exit()
			elif icbc_anal:
				## Build array of date strings for knowing which dates to grab for IC/BC analysis files
				cur_date = beg_date
				icbc_anal_dates = []    # datetime object array
				icbc_anal_dates_str = []   # string array
				icbc_anal_dt = []
				icbc_anal_dt_str = []

				if end_ungrib_exists:
					end_date_icbc = end_date_u
				else:
					end_date_icbc = end_date

				while (cur_date <= end_date_icbc):
					icbc_anal_dates.append(cur_date)
					icbc_anal_dates_str.append(cur_date.strftime(time_fmt_hpss_date_dir))
					cur_date = cur_date + dt.timedelta(days=1)

				cur_dt = beg_datetime
				if icbc_model == 'HRRR' or icbc_model == 'ERA5':
					h_delta = 1
				elif icbc_model == 'GFS' or icbc_model == 'ERA-I':
					h_delta = 6
				else:
					print('ERROR: Unknown icbc_model, cannot set h_delta. Exiting!')
					sys.exit()
				while (cur_dt <= end_datetime):
					icbc_anal_dt.append(cur_dt)
					icbc_anal_dt_str.append(cur_dt.strftime(time_fmt_yyyymmddhh))
					cur_dt = cur_dt + dt.timedelta(hours=h_delta)

				## Download the data from HPSS for entire calendar days (simpler to code, minimal extra time)
				for date_str in icbc_anal_dates_str:
					init_yr	= date_str[0:4]
					## Need to use cget to grab 5 or fewer files at a time to avoid getting locked out by hsi
					if icbc_model == 'RAP':
						print('WARNING: NCAR HPSS will be decommissioned on 1 Oct 2021. A new data source for RAP analyses is required.')
						print('Please edit setup_wrf.py accordingly. Exiting!')
						sys.exit()
#						os.system('hsi cget '+str(icbc_dir)+'/'+date_str+'/*f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+date_str+'/*i0[0-4]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+date_str+'/*i0[5-9]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+date_str+'/*i1[0-4]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+date_str+'/*i1[5-9]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+date_str+'/*i2[0-3]_f000*')
						if ungrib_soil:
							os.system('hsi cget '+str(soil_dir)+'/'+date_str+'/*i0[0-4]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+date_str+'/*i0[5-9]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+date_str+'/*i1[0-4]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+date_str+'/*i1[5-9]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+date_str+'/*i2[0-3]_f000*')

					elif icbc_model == 'HRRR':
						## For downloading from NCAR HPSS prior to 2018 (HPSS end of life on 1 Oct 2021)
						'''
						os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+'/'+date_str+'/*i0[0-4]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+'/'+date_str+'/*i0[5-9]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+'/'+date_str+'/*i1[0-4]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+'/'+date_str+'/*i1[5-9]_f000*')
						os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+'/'+date_str+'/*i2[0-3]_f000*')
						if ungrib_soil:
							os.system('hsi cget '+str(soil_dir)+'/'+init_yr+'/'+date_str+'/*i0[0-4]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+init_yr+'/'+date_str+'/*i0[5-9]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+init_yr+'/'+date_str+'/*i1[0-4]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+init_yr+'/'+date_str+'/*i1[5-9]_f000*')
							os.system('hsi cget '+str(soil_dir)+'/'+init_yr+'/'+date_str+'/*i2[0-3]_f000*')
						'''
						## For downloading from Google Cloud (12 Jul 2018 and later)
						print('WARNING: HRRR section for icbc_anal if branch needs attention')
						fname = 'hrrr.t'+init_hr+'z.wrfnatf00.grib2'
						url = 'https://storage.googleapis.com/high-resolution-rapid-refresh/hrrr.'+date_str+'/'+fname
						if not pathlib.Path(fname).is_file():
							print('   Downloading '+url)
							try:
								wget.download(url)
								print('')
							except:
								wget_error(str(e), now_time_beg)
							
						else:
							print('   File '+fname+' already exists locally. Not downloading again from server.')

						if ungrib_soil:
							fname = 'hrrr.t'+init_hr+'z.wrfprsf00.grib2'
							url = 'https://storage.googleapis.com/high-resolution-rapid-refresh/hrrr.'+date_str+'/'+fname
							if not pathlib.Path(fname).is_file():
								print('   Downloading '+url)
								try:
									wget.download(url)
									print('')
								except Exception as e:
									wget_error(str(e), now_time_beg)
							else:
								print('   File '+fname+' already exists locally. Not downloading again from server.')

					elif icbc_model == 'GFS':
						if gfs_source == 'RDA':
							print('Linking to 0.25-deg GFS output in directory '+str(icbc_dir)+'/'+init_yr+'/'+init_yr+init_mo+init_dy)
							for xx in icbc_anal_dt:
								gfs_file = icbc_dir.joinpath(init_yr, xx.strftime(time_fmt_yyyymmdd), 'gfs.0p25.'+xx.strftime(time_fmt_yyyymmddhh)+'.f000.grib2')
								if gfs_file.is_file():
									os.system('ln -sf '+str(gfs_file))
						elif gfs_source == 'NOMADS':
							print('ERROR: Need to add code to setup_wrf.py to download GFS analyses from NOMADS. Exiting!')
							sys.exit()
						else:
							print('ERROR: gfs_source set to something other than RDA or NOMADS. Please correct this in setup_wrf.py. Exiting!')
							sys.exit()

					elif icbc_model == 'ERA-I':
						print('Linking to ERA-Interim pressure-level output in directory '+str(icbc_dir)+'/'+init_yr+init_mo)
						print('Linking to ERA-Interim surface-level output in directory '+str(surf_dir)+'/'+init_yr+init_mo)
						for xx in icbc_anal_dt:
							this_yr = xx.strftime('%Y')
							this_mo = xx.strftime('%m')
							this_dy = xx.strftime('%d')
							this_hr = xx.strftime('%H')
							erai_pluv_file = icbc_dir.joinpath(this_yr+this_mo, 'ei.oper.an.pl.regn128uv.'+xx.strftime(time_fmt_yyyymmddhh))
							erai_plsc_file	= icbc_dir.joinpath(this_yr+this_mo, 'ei.oper.an.pl.regn128sc.'+xx.strftime(time_fmt_yyyymmddhh))
							erai_surf_file = surf_dir.joinpath(this_yr+this_mo, 'ei.oper.an.sfc.regn128sc.'+xx.strftime(time_fmt_yyyymmddhh))
							if erai_pluv_file.is_file():
								os.system('ln -sf '+str(erai_pluv_file))
							if erai_plsc_file.is_file():
								os.system('ln -sf '+str(erai_plsc_file))
							if erai_surf_file.is_file():
								os.system('ln -sf '+str(erai_surf_file))

					elif icbc_model == 'ERA5':
						print('Linking to ERA5 pressure-level output in directory '+str(icbc_dir)+'/'+init_yr+init_mo)
						print('Linking to ERA5 surface-level output in directory '+str(surf_dir)+'/'+init_yr+init_mo)
						for xx in icbc_anal_dt:
							this_yr = xx.strftime('%Y')
							this_mo = xx.strftime('%m')
							this_dy = xx.strftime('%d')
							this_dir = str(icbc_dir.joinpath(this_yr+this_mo))
							os.system('ln -sf '+str(icbc_dir.joinpath(this_yr+this_mo, 'e5.oper.an.pl.128_*.'+this_yr+this_mo+this_dy+'00_'+this_yr+this_mo+this_dy+'23.grb'))+' .')
							os.system('ln -sf '+str(surf_dir.joinpath(this_yr+this_mo, 'e5.oper.an.sfc.*.grb'))+' .')

			elif icbc_fcst:
				## First parse the icbc forecast cycle string into components
				init_yr	= icbc_fc_cyc[0:4]
				init_mo	= icbc_fc_cyc[5:7]
				init_dy	= icbc_fc_cyc[8:10]
				init_hr	= icbc_fc_cyc[11:13]

				## Build string array of lead hours (restrict to 3-hourly for GFS from NOMADS to limit download/WPS run times, otherwise 1-hourly)
				lead_hrs = ['' for x in range(0, sim_hrs+1, lead_stride)]
				n_lead_hrs = len(lead_hrs)
				for hh in range(n_lead_hrs):
					f_hr = hh*lead_stride + icbc_fc_dt
					if icbc_model == 'HRRR' or icbc_model == 'RAP':
						lead_hrs[hh] = f'{f_hr:02d}'
					elif icbc_model == 'GFS':
						lead_hrs[hh] = f'{f_hr:03d}'
			
				## Need to use cget to grab 5 or fewer files at a time to avoid getting locked out by hsi
				if icbc_model == 'RAP':
					print('WARNING: NCAR HPSS will be decommissioned on 1 Oct 2021. A new data source for RAP forecasts is required.')
					print('Please edit setup_wrf.py accordingly. Exiting!')
					sys.exit()
#					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'*')
					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[0-4]*')
					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[5-9]*')
					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f01[0-4]*')
					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f01[5-9]*')
					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f02[0-4]*')
					if ungrib_soil:
						os.system('hsi cget '+str(soil_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[0-4]*')
						os.system('hsi cget '+str(soil_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[5-9]*')
						os.system('hsi cget '+str(soil_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f01[0-4]*')
						os.system('hsi cget '+str(soil_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f01[5-9]*')
						os.system('hsi cget '+str(soil_dir)+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f02[0-4]*')

				elif icbc_model == 'HRRR':
					## For downloading data from NCAR HPSS prior to 2018 (HPSS end of life on 1 Oct 2021)
					'''
					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[0-4]*')
					os.system('hsi cget '+str(icbc_dir)+'/'+init_yr+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[5-9]*')
					if ungrib_soil:
						os.system('hsi cget '+str(soil_dir)+'/'+init_yr+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[0-4]*')
						os.system('hsi cget '+str(soil_dir)+'/'+init_yr+'/'+init_yr+init_mo+init_dy+'/*i'+init_hr+'_f00[5-9]*')
					'''
					## For downloading data from Google Cloud (data available from 12 Jul 2018 and later)
					for hh in lead_hrs:
						fname = 'hrrr.t'+init_hr+'z.wrfnatf'+hh+'.grib2'
						url = 'https://storage.googleapis.com/high-resolution-rapid-refresh/hrrr.'+init_yr+init_mo+init_dy+'/conus/'+fname
						if not pathlib.Path(fname).is_file():
							print('   Downloading '+url)
							try:
								wget.download(url)
								print('')
							except Exception as e:
								wget_error(str(e), now_time_beg)
						else:
							print('   File '+fname+' already exists locally. Not downloading again from server.')

						if ungrib_soil:
							fname = 'hrrr.t'+init_hr+'z.wrfprsf'+hh+'.grib2'
							url = 'https://storage.googleapis.com/high-resolution-rapid-refresh/hrrr.'+init_yr+init_mo+init_dy+'/conus/'+fname
							if not pathlib.Path(fname).is_file():
								print('   Downloading '+url)
								try:
									wget.download(url)
									print('')
								except Exception as e:
									wget_error(str(e), now_time_beg)
							else:
								print('   File '+fname+' already exists locally. Not downloading again from server.')

				elif icbc_model == 'GFS':
					if gfs_source == 'RDA':
						print('Linking to 0.25-deg GFS output in directory '+str(icbc_dir)+'/'+init_yr+'/'+init_yr+init_mo+init_dy)
					elif gfs_source == 'NOMADS':
						print('Downloading 0.25-deg GFS output from directory '+str(icbc_dir)+'/'+init_yr+'/'+init_yr+init_mo+init_dy)

					for hh in lead_hrs:
						if gfs_source == 'RDA':
							## If linking to GFS on GLADE RDA:
							gfs_file = icbc_dir.joinpath(init_yr, init_yr+init_mo+init_dy, 'gfs.0p25.'+init_yr+init_mo+init_dy+init_hr+'.f'+hh+'.grib2')
							if gfs_file.is_file():
								os.system('ln -sf '+str(gfs_file))
						elif gfs_source == 'NOMADS':
							## If downloading GFS from NCEP NOMADS:
							fname = 'gfs.t'+init_hr+'z.pgrb2.0p25.f'+hh
							url = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.'+init_yr+init_mo+init_dy+'/'+init_hr+'/'+fname
							if not pathlib.Path(fname).is_file():
								print('   Downloading '+url)
								try:
									wget.download(url)
									print('')
								except Exception as e:
									wget_error(str(e), now_time_beg)
							else:
								print('   File '+fname+' already exists locally. Not downloading again from server.')
						else:
							print('ERROR: gfs_source set to something other than RDA or NOMADS. Please correct this in setup_wrf.py. Exiting!')
							sys.exit()

			else:
				print('   ERROR: Both icbc_anal and icbc_fcst = False. Set one to True and re-run this script. Exiting!')
				sys.exit()

			## Unzip the files
			if icbc_model == 'RAP':
				os.system('gunzip *grb2.gz')
			os.chdir(str(run_dir))

		if (icbc_anal or icbc_fcst) and not get_icbc:
			print('WARNING: get_icbc = False, so even though icbc_anal or icbc_fcst are True, no IC/BC grib files will be downloaded.')
	

		## Modify the namelist for running ungrib and/or avg_tsfc with the selected start/end dates
		if run_ungrib and run_avg_tsfc:
			if end_ungrib_exists:
				end_str_ungrib = end_ungrib
			else:
				end_str_ungrib = end_dt

			with open('namelist.wps', 'r') as input_file, open('namelist.wps.ungrib', 'w') as output_file:
				for line in input_file:
					if line.strip()[0:8] == 'end_date':
						output_file.write(" end_date   = '"+end_str_ungrib+":00:00',\n")
					else:
						output_file.write(line)
		
		## Run ungrib
		if run_ungrib:
			ungrib_dir.mkdir(parents=True, exist_ok=True)

			if ungrib_soil:
				print('Running ungrib for soil and/or surface variables . . .')
				if icbc_model == 'HRRR':
					os.system('./link_grib.csh '+str(grib_dir)+'/*wrfprs*.gr*b2')
				else:
					print('ERROR: Need to set link_grib directives for this icbc_model in ungrib_soil. Exiting!')
					sys.exit()

				if pathlib.Path('Vtable').is_file():
					pathlib.Path('Vtable').unlink()
				pathlib.Path('Vtable').symlink_to(wps_dir.joinpath('ungrib','Variable_Tables',soil_vtable))
				if pathlib.Path('namelist.wps').is_file():
					pathlib.Path('namelist.wps').unlink()
				pathlib.Path('namelist.wps').symlink_to(pathlib.Path(soil_nml))
				os.system('./ungrib.exe >& ungrib.out.soil')
				if 'Successful completion' in open('ungrib.log').read():
					print('   SUCCESS! ungrib completed successfully.')
					shutil.move('ungrib.log','ungrib.log.soil')
				else:
					print('   ERROR: ungrib.exe failed.')
					print('   Consult '+str(run_dir)+'/ungrib.log.soil and ungrib.out.soil for potential error messages.')
					print('   Exiting !')
					shutil.move('ungrib.log','ungrib.log.soil')
					sys.exit()

			print('Running ungrib for IC/BC data . . .')
			for x in pathlib.Path('.').glob('GRIBFILE.*'):
				pathlib.Path(x).unlink()
			if icbc_model == 'HRRR':
				os.system('./link_grib.csh '+str(grib_dir)+'/*wrfnat*.gr*b2')
			elif icbc_model == 'GFS':
				os.system('./link_grib.csh '+str(grib_dir)+'/gfs.*')
			elif icbc_model == 'ERA-I':
				os.system('./link_grib.csh '+str(grib_dir)+'/ei.oper.an.*')
			elif icbc_model == 'ERA5':
				os.system('./link_grib.csh '+str(grib_dir)+'/e5.oper.an.*')
			else:
				print('ERROR: Need to set link_grib directives for this icbc_model. Exiting!')
				sys.exit()

			if pathlib.Path('Vtable').is_file():
				pathlib.Path('Vtable').unlink()
			pathlib.Path('Vtable').symlink_to(wps_dir.joinpath('ungrib','Variable_Tables',icbc_vtable))
			if run_avg_tsfc:
				if pathlib.Path('namelist.wps').is_file():
					pathlib.Path('namelist.wps').unlink()
				pathlib.Path('namelist.wps').symlink_to(pathlib.Path('namelist.wps.ungrib'))
			else:
				if 'namelist.wps' != icbc_nml:	# also could perhaps use the condition 'if ungrib_soil'
					## If icbc_nml is already namelist.wps, then we don't want to link namelist.wps to itself, because it's already linked properly
					if pathlib.Path('namelist.wps').is_file():
						pathlib.Path('namelist.wps').unlink()
					pathlib.Path('namelist.wps').symlink_to(pathlib.Path(icbc_nml))
				os.system('./ungrib.exe >& ungrib.out')
			if 'Successful completion' in open('ungrib.log').read():
				print('   SUCCESS! ungrib completed successfully.')
			else:
				print('   ERROR: ungrib.exe failed.')
				print('   Consult '+str(run_dir)+'/ungrib.log and ungrib.out for potential error messages.')
				print('   Exiting!')
				sys.exit()

			if run_avg_tsfc:
				if pathlib.Path('namelist.wps').is_file():
					pathlib.Path('namelist.wps').unlink()
				pathlib.Path('namelist.wps').symlink_to(pathlib.Path(arc_dir.joinpath('config','namelist.wps')))

		## Run avg_tsfc
		if run_avg_tsfc:
			ungrib_dir.mkdir(parents=True, exist_ok=True)
			if pathlib.Path('namelist.wps').is_file():
				pathlib.Path('namelist.wps').unlink()
			pathlib.Path('namelist.wps').symlink_to(pathlib.Path('namelist.wps.ungrib'))
			print('Running avg_tsfc . . .')
			os.system('./avg_tsfc.exe >& avg_tsfc.out')
			if 'Successful completion' in open('avg_tsfc.out').read():
				print('   SUCCESS! avg_tsfc completed successfully.')
				shutil.move('TAVGSFC',str(ungrib_dir))
				pathlib.Path('TAVGSFC').symlink_to(ungrib_dir.joinpath('TAVGSFC'))
				if pathlib.Path('namelist.wps').is_file():
					pathlib.Path('namelist.wps').unlink()
				pathlib.Path('namelist.wps').symlink_to(arc_dir.joinpath('config','namelist.wps'))
			else:
				print('   ERROR: avg_tsfc.exe failed.')
				print('   Consult '+str(run_dir)+'/avg_tsfc.out for potential error messages.')
				print('   Exiting!')
				sys.exit()

		## Run metgrid
		if run_metgrid:
			met_dir.mkdir(parents=True, exist_ok=True)
			if not pathlib.Path('namelist.wps').is_file():
				if pathlib.Path('namelist.wps.hybr').is_file():
					pathlib.Path('namelist.wps').symlink_to(pathlib.Path('namelist.wps.hybr'))
				else:
					print('ERROR: Cannot find a file to link to from namelist.wps to run metgrid. Exiting!')
					sys.exit()
			print('Running metgrid . . .')
			os.system('./metgrid.exe >& metgrid.out')
			if 'Successful completion' in open('metgrid.log').read():
				print('   SUCCESS! metgrid completed successfully.')
				if 'WARNING' in open('metgrid.out').read():
					print('      WARNING messages detected in '+str(run_dir)+'/metgrid.out. Check the file to investigate.')
			else:
				print('   ERROR: metgrid.exe failed.')
				print('   Consult '+str(run_dir)+'/metgrid.log and metgrid.out for potential error messages.')
				print('   Exiting!')
				sys.exit()

		## Submit real
		if sub_real:
			print('Submitting real to the queue . . .')
			if pathlib.Path('rsl.out.0000').is_file():
				for x in pathlib.Path('.').glob('rsl.*'):
					pathlib.Path(x).unlink()
			for x in pathlib.Path('.').glob('met_em*'):
				pathlib.Path(x).unlink()
			for x in pathlib.Path(met_dir).glob('met_em*'):
				pathlib.Path(x.name).symlink_to(x)
			os.system('qsub pbs_submit_real.bash')

			## Monitor the progress of real
			status = False
			while not status:
				if not os.path.exists('rsl.out.0000'):
					time.sleep(5)
				else:
					print('real is now running . . .')
					status = True
			status = False
			while not status:
				if 'real_em: SUCCESS COMPLETE REAL_EM INIT' in open('rsl.out.0000').read():
					print('   SUCCESS! real completed successfully.')
					for x in pathlib.Path('.').glob('rsl.*'):
						shutil.move(str(x), run_dir.joinpath('rsl_real', x))
					time.sleep(1)	# Brief pause to let the file system gather itself (found a problem here with cron)
					status = True
				else:
					if 'FATAL' in open('rsl.out.0000').read():
						print('   ERROR: real.exe failed.')
						print('   Consult '+str(run_dir)+'/rsl.out.0000 for potential error messages.')
						print('   Exiting!')
						sys.exit()
					time.sleep(5)

		## Copy aerosol forecast data from another server
		if aero_scp_data:
			os.chdir(aero_dir)
			aero_scp_path_dir	= aero_scp_path.joinpath(beg_dt_yyyymmddhh)

			## Don't bother scp'ing files if they already exist locally
			if aero_dir.joinpath(beg_dt_yyyymmddhh).is_dir():
				print(aero_model+' files already exist locally. Not copying files from remote server.')
				print('Local dir: '+str(aero_dir.joinpath(beg_dt_yyyymmddhh)))
			else:
				print('Copying files from '+aero_scp_machine+':'+str(aero_scp_path_dir)+' to '+str(aero_dir))
				os.system('scp -r '+aero_scp_machine+':'+str(aero_scp_path_dir)+' .')

			## If the scp wasn't successful, then look for a prior time
			if not aero_dir.joinpath(beg_dt_yyyymmddhh).is_dir():
				scp_success = False
				try_this_dt = beg_datetime
				if aero_model == 'GEOS-5':
					## GEOS-5 forecasts go out to 240 h for 00z cycles, 120 h for 12z cycles, and 30 h for 06z & 18z cycles
					for xx in range(1, 241):
						try_this_dt = try_this_dt - dt.timedelta(hours=xx)
						try_this_dt_str = try_this_dt.strftime(time_fmt_yyyymmddhh)
						if not (try_this_dt.hour == 0 or try_this_dt.hour == 6 or try_this_dt.hour == 12 or try_this_dt.hour == 18):
							continue
						elif (try_this_dt.hour == 6 or try_this_dt.hour == 18) and (sim_hrs + xx) > 30:	# The 06z & 18z GEOS-5 cycles run out to 30 h
							continue
						elif try_this_dt.hour == 12 and (sim_hrs + xx) > 120:		# The 12z GEOS-5 cycle runs out to 120 h
							continue
						elif try_this_dt.hour == 0 and (sim_hrs + xx) > 240:		# The 00z GEOS-5 cycle runs out to 240 h
							print('ERROR: No available GEOS-5 forecasts can be copied from '+aero_scp_machine+' that span the WRF simulation time. Exiting!')
							sys.exit()
						else:
							aero_scp_path_dir = aero_scp_path.joinpath(try_this_dt_str)
							## Don't bother scp'ing files if they already exist locally
							if aero_dir.joinpath(try_this_dt_str).is_dir():
								print('Found GEOS-5 files that already exist locally. Not copying files from remote server.')
								print('Local dir: '+str(aero_dir.joinpath(try_this_dt_str)))
								break
							else:
								print('Copying files from '+aero_scp_machine+':'+str(aero_scp_path_dir)+' to '+str(aero_dir))
								os.system('scp -r '+aero_scp_machine+':'+str(aero_scp_path_dir)+' .')

								## If this scp also wasn't successful, then look further back for another prior time
								if aero_dir.joinpath(try_this_dt_str).is_dir():
									break		# success!
								else:
									continue
				elif aero_model == 'CAMS':
					## CAMS forecasts are issued at 00 and 12 UTC, with hourly output out to 120 h
					## NOTE: This section of code has not yet been fully tested
					for xx in range(1,121):
						try_this_dt = try_this_dt - dt.timedelta(hours=xx)
						try_this_dt_str = try_this_dt.strftime(time_fmt_yyyymmddhh)
						if not (try_this_dt.hour == 0 or try_this_dt.hour == 12):
							continue
						else:
							aero_scp_path_dir = aero_scp_path.joinpath(try_this_dt_str)
							## Don't bother scp'ing files if they already exist locally
							if aero_dir.joinpath(try_this_dt_str).is_dir():
								print('Found CAMS files that already exist locally. Not copying files from remote server.')
								print('Local dir: '+str(aero_dir.joinpath(try_this_dt_str)+' .'))
								break
							else:
								print('Copying files from '+aero_scp_machine+':'+str(aero_scp_path_dir)+' to '+str(aero_dir))
								os.system('scp -r '+aero_scp_machine+':'+str(aero_scp_path_dir)+' .')

								## If this scp also wasn't successful, then look further back for another prior time
								if aero_dir.joinpath(try_this_dt_str).is_dir():
									break		# success!
								else:
									continue

		## Impose aerosols (AOD) onto WRF
		## Currently this is only supported from GEOS-5
		## NOTE: This is a re-creation of the functionality in: /glade/work/jimenez/devel/wrfauto/io/impose_aerosols.s
		if impose_aerosol:
			os.chdir(run_dir)

			## Clean up previous runs
			if pathlib.Path('aerinput_d01').is_file():
				for x in pathlib.Path('.').glob('aerinput_d0*'):
					pathlib.Path(x).unlink()

			## Calculate number of files to process
			n_files_to_process = sim_hrs // aero_freq_h + 1

			## Check for aerosol files to determine if imposing aerosols is possible
			## NOTE: This portion is a re-creation of the functionality in: /glade/work/jimenez/devel/wrfauto/io/check_aerosols_geos.s
			is_sim_possible = False
			aero_init_datetime = beg_datetime
			for hh in range(0, (240-sim_hrs+1)):
				aero_init_datetime = beg_datetime - dt.timedelta(hours=hh)
				aero_init_dt_str_dir = aero_init_datetime.strftime(time_fmt_yyyymmddhh)
				aero_init_dt_str_fil = aero_init_datetime.strftime(time_fmt_yyyymmdd_hh)
				aero_dir_mostrecent = aero_dir.joinpath(aero_init_dt_str_dir)

				## Only bother looking for the directory and files if it's possible it exists and can meet our forecast requirements
				if not (aero_init_datetime.hour == 0 or aero_init_datetime.hour == 6 or aero_init_datetime.hour == 12 or aero_init_datetime.hour == 18):
					continue
				elif (aero_init_datetime.hour == 6 or aero_init_datetime.hour == 18) and (sim_hrs+hh > 30):	# 06z and 18z GEOS-5 runs only go to 30 h
					continue
				elif aero_init_datetime.hour == 12 and sim_hrs+hh > 120:	# 12z GEOS-5 runs go to 120 h
					continue
				elif sim_hrs+hh > 240:	# the longest GEOS-5 runs are the 00z cycles that go to 240 h, so if we get here, it's not possible to impose aerosols
					break

				## Now see if there are enough files to successfully impose aerosols
				if aero_dir_mostrecent.is_dir():
					file_list = sorted(aero_dir_mostrecent.glob('GEOS.fp.fcst.inst1_2d_hwl_Nx.'+aero_init_dt_str_fil+'+*.V01.nc4'))
					try:
						ind_beg	= file_list.index( \
										aero_dir_mostrecent.joinpath('GEOS.fp.fcst.inst1_2d_hwl_Nx.'+aero_init_dt_str_fil+'+'+beg_dt_yyyymmdd_hhmm+'.V01.nc4'))
					except ValueError:
						ind_beg	= None
					try:
						ind_end	= file_list.index( \
										aero_dir_mostrecent.joinpath('GEOS.fp.fcst.inst1_2d_hwl_Nx.'+aero_init_dt_str_fil+'+'+end_dt_yyyymmdd_hhmm+'.V01.nc4'))
					except ValueError:
						ind_end	= None
			
					## We need GEOS-5 files every hour through the duration of the WRF simulation in order to impose aerosols
					if ind_beg != None and ind_end != None and ind_end-ind_beg+1 == n_files_to_process:
						is_sim_possible = True
						break
					else:
						continue		# keep looking
				else:
					continue		# keep looking
	
			if is_sim_possible:
				## Convert GEOS-5 aerosol files to WRF grid and expected format
				## Link to the aerosol files
				print('Linking to aerosol files in directory '+str(aero_dir_mostrecent))
				for aero_fname in file_list[ind_beg:ind_end+1]:
					if pathlib.Path(aero_fname.name).is_symlink():
						pathlib.Path(aero_fname.name).unlink()
					pathlib.Path(aero_fname.name).symlink_to(aero_fname)

				## Loop over each domain
				for dd in range(1, num_dom+1):
					## First write text file to be piped as input to geos2wrf.exe
					with open('aero.txt', 'w') as outfile:
						outfile.write('%s' % 'wrfinput_d0'+str(dd)+'\n')
						outfile.write('%s' % str(n_files_to_process)+'\n')
						for aero_fname in file_list[ind_beg:ind_end+1]:
							outfile.write('%s\n' % aero_fname.name)
						outfile.write('%s' % 'aerinput_d0'+str(dd)+'\n')
					print('Running '+str(aero_exe)+' for domain '+str(dd)+' files...')
					os.system(str(aero_exe)+' < aero.txt')

			else:
				print('ERROR: No aerosol data found that will cover the entire WRF simulation period.')
				print('       Possible solutions:')
				print('       1. Check the aerosol data source')
				print('       2. Change the simulation start date and/or simulation length to accomodate available aerosol data')
				print('       3. Change the default WRF namelist to not expect aerinput file(s), and re-run this script with impose_aerosol = False')
				print('Exiting!')
				sys.exit()

		## Submit wrf
		if sub_wrf:
			os.chdir(run_dir)
			if not is_skebs_ens:
				print('Submitting wrf to the queue . . .')
				os.system('qsub pbs_submit_wrf.bash')

				## Monitor the progress of wrf
				if rt_flag:
					status = False
					while not status:
						if not os.path.exists('rsl.out.0000'):
							time.sleep(10)
						else:
							print('wrf is now running . . .')
							status = True
					status = False
					while not status:
						if 'wrf: SUCCESS COMPLETE WRF' in open('rsl.out.0000').read():
							print('   SUCCESS! WRF completed successfully.')
							for x in pathlib.Path('.').glob('rsl.*'):
								shutil.move(str(x), run_dir.joinpath('rsl_wrf', x))
							time.sleep(1)	# Brief pause to let the file system gather itself (found a problem here with cron)
							status = True
						else:
							if 'FATAL' in open('rsl.out.0000').read():
								print('   ERROR: wrf.exe failed.')
								print('   Consult '+str(run_dir)+'/rsl.out.0000 for potential error messages.')
								print('   Exiting!')
								sys.exit()
							time.sleep(5)
			else:
				for member in ens_mem_arr_str:
					os.chdir(str(run_dir)+'/mem'+member)
					print('Submitting wrf to the queue for SKEBS ensemble member '+member+' . . .')

					## All SKEBS ensemble members use the same wrfinput & wrfbdy files
					for x in pathlib.Path('.').glob('wrfinput*'):
						pathlib.Path(x).unlink()
					for x in pathlib.Path('..').glob('wrfinput*'):
						pathlib.Path(x.name).symlink_to(x)
					if pathlib.Path('wrfbdy_d01').is_file():
						pathlib.Path('wrfbdy_d01').unlink()
					pathlib.Path('wrfbdy_d01').symlink_to(pathlib.Path('../wrfbdy_d01'))

					os.system('qsub pbs_submit_wrf.bash')
					os.chdir(str(run_dir))

		## Archive wrf files
		if arc_wrf:
			os.chdir(run_dir)
			print('Moving wrfinput/wrfbdy files to '+str(arc_dir)+'/config')
			for x in pathlib.Path('.').glob('wrfinput*'):
				shutil.move(str(x),str(arc_dir.joinpath('config')))
			for x in pathlib.Path('.').glob('wrfbdy*'):
				shutil.move(str(x),str(arc_dir.joinpath('config')))
			if arc_gzip:
				os.system('gzip '+str(arc_dir)+'/config/wrf*')

			if not is_skebs_ens:
				if arc_wrfout:
					print('Moving wrfout files to '+str(arc_dir)+'/wrfout')
					for x in pathlib.Path('.').glob('wrfout*'):
						shutil.move(str(x), str(arc_dir.joinpath('wrfout')))
					if arc_gzip:
						os.system('gzip '+str(arc_dir)+'/wrfout/*')
				if arc_uppout:
					print('Moving UPP grib2 output files to '+str(arc_dir)+'/uppout')
					os.chdir(upp_run_dir)
					for x in pathlib.Path('.').glob('WRFPRS*.grib2'):
						shutil.move(str(x),str(arc_dir.joinpath('uppout')))

				if write_ts:
					print('Moving time series files to '+str(arc_dir)+'/ts')
					for x in pathlib.Path('.').glob('*.TS'):
						shutil.move(str(x), str(arc_dir.joinpath('ts')))
					for x in pathlib.Path('.').glob('*.PH'):
						shutil.move(str(x), str(arc_dir.joinpath('ts')))
					for x in pathlib.Path('.').glob('*.QV'):
						shutil.move(str(x), str(arc_dir.joinpath('ts')))
					for x in pathlib.Path('.').glob('*.TH'):
						shutil.move(str(x), str(arc_dir.joinpath('ts')))
					for x in pathlib.Path('.').glob('*.UU'):
						shutil.move(str(x), str(arc_dir.joinpath('ts')))
					for x in pathlib.Path('.').glob('*.VV'):
						shutil.move(str(x), str(arc_dir.joinpath('ts')))
					if arc_gzip:
						os.system('gzip '+str(arc_dir)+'/ts/*')
			else:
				for member in ens_mem_arr_str:
					os.chdir(str(run_dir)+'/mem'+member)
					if arc_wrfout:
						print('Moving wrfout files to '+str(arc_dir)+'/mem'+member+'/wrfout')
						for x in pathlib.Path('.').glob('wrfout*'):
							shutil.move(str(x), str(arc_dir.joinpath('mem'+member,'wrfout')))
						if arc_gzip:
							os.system('gzip '+str(arc_dir)+'/mem'+member+'/wrfout/*')

if __name__ == '__main__':
	now_time_beg = dt.datetime.utcnow()
	init_dt_first, init_dt_last, init_stride_h, sim_hrs, icbc_fc_dt, rt_flag, arc_flag = parse_args()
	main(init_dt_first, init_dt_last, init_stride_h, sim_hrs, icbc_fc_dt, rt_flag, arc_flag, now_time_beg)
	now_time_end = dt.datetime.utcnow()
	run_time_tot = now_time_end - now_time_beg
	now_time_beg_str = now_time_beg.strftime('%Y-%m-%d %H:%M:%S')
	now_time_end_str = now_time_end.strftime('%Y-%m-%d %H:%M:%S')
	print('\nScript completed successfully.')
	print('   Beg time: '+now_time_beg_str)
	print('   End time: '+now_time_end_str)
	print('   Run time: '+str(run_time_tot)+'\n')
