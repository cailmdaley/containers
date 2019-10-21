from spt3g import core, coordinateutils, mapmaker
import pickle

m = list(core.G3File('/home/nlharr/tmp/maps500d/sptpol_season.g3'))[0]['T']
m *= 0
mapmaker.pointsourceutils.make_point_source_map(
    m, '/home/nlharr/tmp/maps500d/ptsrc_config_ra0hdec-57p5_both_50mJy.txt')
pickle.dump(m, open('point_source_mask.pkl', 'w'))

