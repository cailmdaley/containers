from spt3g import core, sptpol

reader = sptpol.SPTDataReader(master_configfile='sptpol_master_config', start_date='130406 05:24:18', stop_date='130406 05:28:21')
reader.readData(verbose=True, zero_processing=False, remove_unknown_bolos=True, timestream_units = 'amps')

pipe = core.G3Pipeline()
pipe.Add(reader)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer('scans.g3'))

pipe.Run()

