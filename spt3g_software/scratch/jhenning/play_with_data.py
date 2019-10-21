from spt3g import core, gcp, std_processing, dfmux, coordinateutils, sptpol


p = core.G3Pipeline()
p.Add(core.G3Reader, filename='/home/jhenning/data/RCW38/882463/0002.g3')

pointing_info = []
def grabTrackerPointing(fr):
    if 'TrackerPointing' in fr:
        print 'We have printing!'
        pointing_info.append(fr['TrackerPointing'])

p.Add(grabTrackerPointing)
p.Run()

