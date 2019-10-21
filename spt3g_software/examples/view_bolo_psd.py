from spt3g import core, mapmaker, sptpol
from spt3g.todfilter import ViewTimestreamPSD

bolos_to_disp = ['C1.A1.1.X','C1.A1.1.Y', 'C1.A1.10.X', 'C1.A1.10.Y',
                 'C1.A1.11.X', 'C1.A1.11.Y', 'C1.A1.12.X', 'C1.A1.12.Y']

pipe = core.G3Pipeline()
idf_fn = mapmaker.get_test_files_path() + 'ra23h30dec-55_idf_20130430_003228_150ghz.h5'
pipe.Add(sptpol.DirectIdfReader, filename = idf_fn)

# Low frequency cutoff to eliminate 1/f noise
low_f = 1.8*core.G3Units.Hz

# With bolo_order present, plots each listed bolo in a separate subplot.
# pipe.Add(ViewTimestreamPSD, bolo_order=bolos_to_disp, low_f=low_f)
pipe.Add(ViewTimestreamPSD, low_f=low_f)

pipe.Run()
