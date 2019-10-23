# dfmux_process.py
#
# Script for doing processing of dfmux data. By 'processing', we merely
# mean collating the raw data stream into scans and adding some additional
# data to those scan frames. At the moment, these are mainly pointing 
# information and scan flags from the SPTpol arcfiles. All of the work is
# actually done by functions within the 'dfmuxtools' python module.
#
# This is a dumbed-down version of an previous attempt at a processing
# script. Since we are only usually processing one GCP observation at a
# time with the IceBoard in our current testing regime, we can just edit
# the observation boundaries by hand in this script.
#
# Adam Anderson
# adama@fnal.gov

import dfmuxtools
import datetime as dt

arcfile_path = '/data/sptdaq/iceboard_testing/arcfiles/'
rawdata_path = '/data/sptdaq/iceboard_testing/dfmux/raw'
procdata_path = '/data/sptdaq/iceboard_testing/dfmux/processed/'
log_path = '/data/sptdaq/iceboard_testing/gcplog/2015/'
time_offset = 29907677.  # see docstring for dfmuxtools.process_observation()
scan_type = 'GCP'        # see docstring for dfmuxtools.process_observation()


#obs_times = [[dt.datetime(2015, 12, 13, 4, 35, 0), dt.datetime(2015, 12, 13, 5, 25, 0)]]
#raw_files = [['/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_043447_00.g3', 
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_043447_01.g3']]

#obs_times = [[dt.datetime(2015, 12, 13, 6, 30, 0), dt.datetime(2015, 12, 13, 7, 20, 0)]]
#raw_files = [['/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_061542_00.g3', 
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_061542_01.g3',
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_061542_02.g3',
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_061542_03.g3',
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_061542_04.g3']]

# attempt to concatenate two separate runs of the DAQ with identical settings
#obs_times = [[dt.datetime(2015, 12, 13, 8, 20, 0), dt.datetime(2015, 12, 13, 9, 52, 0)]]
#raw_files = [['/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_081658_00.g3',
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_081658_01.g3',
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_091229_00.g3',
#              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_091229_01.g3']]

#obs_times = [[dt.datetime(2015, 12, 13, 9, 57, 35), dt.datetime(2015, 12, 13, 12, 0, 41)]]
obs_times = [[dt.datetime(2015, 12, 13, 12, 14, 2), dt.datetime(2015, 12, 13, 14, 23, 9)]]
raw_files = [['/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_00.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_01.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_02.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_03.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_04.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_05.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_06.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_07.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_08.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_09.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_10.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_11.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_12.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_13.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_14.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_15.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_16.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_17.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_18.g3',
              '/data/sptdaq/iceboard_testing/dfmux/raw/dfmux_20151213_095543_19.g3']]

for jobs,obs in enumerate(obs_times):
    dfmuxtools.process_observation(obs[0], obs[1], time_offset, raw_files[jobs], arcfile_path, procdata_path, scan_type)
