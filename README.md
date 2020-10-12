# setup_wrf_epri
Python script to set up and run all the steps of WPS and WRF on NCAR's Cheyenne HPC, with an optional Python wrapper script to enable easier execution from crontab.

setup_wrf.py:
=============
This script links to (and downloads, if necessary) all the files needed to run WPS/WRF.
Each program in the WPS/WRF workflow can be optionally executed if its option is set to True.
WRF output files can also be optionally moved to an archival directory (arc_dir).

usage:

      setup_wrf.py [-h] [-l INIT_DT_LAST] [-i INIT_STRIDE_H] [-s SIM_LENGTH_H] [-f ICBC_FC_DT] [-r] [-a] init_dt_first

positional arguments:

init_dt_first

      beginning date/time of first WRF simulation [YYYYMMDD_HH]

optional arguments:

-h, --help

      show this help message and exit

-l INIT_DT_LAST, --init_dt_last INIT_DT_LAST

      beginning date/time of last WRF simulation [YYYYMMDD_HH] (default: same as init_dt_first)

-i INIT_STRIDE_H, --init_stride_h INIT_STRIDE_H

      integer number of hours between forecast cycles (default: 24)

-s SIM_LENGTH_H, --sim_length_h SIM_LENGTH_H

      integer number of hours for WRF simulation (default: 48)

-f ICBC_FC_DT, --icbc_fc_dt ICBC_FC_DT

      integer number of hours prior to WRF init time for IC/LBC model cycle (default: 0)

-r, --realtime

      flag when running in real-time to keep this script running until WRF is done

-a, --archive

      flag to archive wrfout, wrfinput, wrfbdy, and namelist files off of scratch


setup_wrf_wrapper.py
====================
This script is intended to be called by a crontab job with 3 arguments: WRF start hour, WRF simulation length, and IC/LBC delta hours.
This script gets the current UTC date to pair with the specified WRF start hour to provide the necessary arguments for calling setup_wrf.py.

usage:

      setup_wrf_wrapper.py [-h] [-s SIM_LENGTH_H] [-f ICBC_FC_DT] init_hr

positional arguments:

init_hr

      two-digit hour of WRF start time (e.g., 06)

optional arguments:

  -h, --help

      show this help message and exit

  -s SIM_LENGTH_H, --sim_length_h SIM_LENGTH_H

      integer number of hours for WRF simulation (default: 6)

  -f ICBC_FC_DT, --icbc_fc_dt ICBC_FC_DT

      integer number of hours (default: 0
