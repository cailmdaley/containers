from spt3g import core, mapmaker, sptpol
from spt3g.util import G3Viewer

idf_fn = mapmaker.get_test_files_path() + 'ra23h30dec-55_idf_20130430_003228_150ghz.h5'
reader = sptpol.DirectIdfReader(filename=idf_fn)

frames = reader(None)

v = G3Viewer(frames, current_year=1992)
