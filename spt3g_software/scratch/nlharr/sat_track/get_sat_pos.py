from spt3g import core, util
import numpy as np
import astropy.coordinates, astropy.units, astropy.time
import ephem, glob
import os.path

GU = util.genericutils

'''
(-16100450, -186.34780092592592, 62, 1680, -0.2389427596623751, -0.9487288730917052, 'Sq5SBpol28Ch7', '2016-06-28T16:16:15.248059850')
(-16100450, -186.34780092592592, 62, 1737, -0.23520737356087393, -0.9449650042173952, 'Sq1SBpol18Ch11', '2016-06-28T16:16:15.546904730')
(-16100450, -186.34780092592592, 62, 1781, -0.23257740261985538, -0.9422059703943717, 'Sq3SBpol17Ch5', '2016-06-28T16:16:15.777592010')
(-16100450, -186.34780092592592, 62, 1805, -0.23107968124696696, -0.9406465834613588, 'Sq6SBpol17Ch2', '2016-06-28T16:16:15.903421440')
'''



#generate list of events, need times, ras, decs
#

#curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=NLHARR@berkeley.edu&password=asdfasdfasdfasdf'
#curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2016-06-27--2016-06-29/format/3le' > 3le1606
#curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org//ajaxauth/logout


def t_to_date_str(t_str):
    t_base = core.G3Time(t_str)
    t_low = core.G3Time()
    t_low.mjd = t_base.mjd - 1
    t_high = core.G3Time()
    t_high.mjd = t_base.mjd + 1
    
    ds = [ t_low.isoformat().split('T')[0],
           t_high.isoformat().split('T')[0] ]
    date_str = '%s--%s'%(ds[0], ds[1])
    return date_str

def from_time_get_space_track_request(t_str):
    date_str = t_to_date_str(t_str)
    return """curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/%s/format/3le' > 3le_%s.tles""" %(date_str, date_str), '3le_%s.tles'%date_str


def generate_satellite_list(
        fn, root_path = '/home/nlharr/spt3g_software/scratch/nlharr/sat_track/'):
    #print('reading ', fn)
    f = open(root_path + fn)
    tles = list(f)
    f.close()
    #print('done reading')

    out_sats = []

    if len(tles)%3 != 0:
        import pdb; pdb.set_trace()
    #print("Generating sat lst", len(tles))
    for i in range(0,len(tles), 3):
        out_sats.append(ephem.readtle(tles[i], tles[i+1], tles[i+2]))
    return out_sats

def get_observed_ra_dec(sats, t_string, ev_ra, ev_dec, do_parkes = False, 
                        offset_ra = 0,
                        dist = 0.03*1.0):

    ev_ra = (ev_ra + offset_ra)%(2*np.pi)

    if do_parkes:
        assert(False)
        spt = ephem.Observer()
        spt.lat= '-32.9983'
        spt.lon= '148.2636'
        spt.elevation= 414.8
        spt.date = ephem.Date(t_string.replace('T', ' '))
    else:
        spt = ephem.Observer()
        spt.lat= '-89.991066'
        spt.lon= '-44.65'
        spt.elevation= 2835.0
        spt.date = ephem.Date(t_string.replace('T', ' '))

    out_vals = []

    for i, s in enumerate(sats):
        try:
            s.compute(spt)
            ra = float(s.ra)
            dec = float(s.dec)
            if ((abs(ra - ev_ra)/np.cos(ev_dec) < dist) and
                (abs(dec - ev_dec) < dist)):
                print(s)
                print(i, float(s.ra), float(s.dec))
                out_vals.append( (s.name, float(s.ra), float(s.dec)) )
        except:
            continue
    return out_vals


#high_sig_one.txt  multiple_coincident.txt  two_coincident.txt
#fn = 'evs_lst/high_sig_one.txt'
#fn = 'evs_lst/multiple_coincident.txt'
#fn = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/other_plots/ldfs_good_counted_ldfs_S4p0W4p0H1p0Nv0p0Xv2p0E-04P5p0E-03Ng0Gr2Co7p0Fs1_trip.txt'

#fn = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/other_plots/ldfs_phinal_ldfs_S3p0_2p0W4p0H0p5Nv0p0Xv1p5E-04P1p5E-03Ng0Gr2Co6p5Fs1_trip.txt'

#fn = 'evs_lst/two_coincident.txt'

#fn = 'evs_lst/sig_9_trips.txt'
#fn = 'evs_lst/sig7_trips.txt'
#fn = 'evs_lst/offset_req.txt'

#fn = 'evs_lst/parkes_frb.txt'
if False:
    fn = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/other_plots/ldfs_phinal_ldfs_S3p0_2p0W4p0H0p5Nv0p0Xv1p5E-04P1p3E-03Ng0Gr2Co6p5Fs1_trip.txt'
    f = open(fn)
    info_triplets_inner = map(eval, list(f))
    f.close()
else:
    #fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/other_plots/ldfs_phinal_NEW_GLTCH_G99_p75_ldfs_*')
    #fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/other_plots/ldfs_phinal_NEW_ALL_SQ_GLTCH_G99_*')
    fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/other_plots/ldfs_phinal_phinal_v2_G98_*')
    info_triplets = []
    for fn in fns:
        f = open(fn)
        info_triplets_inner = map(eval, list(f))
        f.close()
        info_triplets += info_triplets_inner
    info_triplets = GU.uniquify_list(map(tuple, info_triplets))


def check_lst(l, bset = ['0 GRACE 2', '0 GRACE 1'],
              #bset = ['0 SWARM B', '0 SWARM A', '0 SWARM C'], 
              do_all = False):
    if do_all:
        return len(l) == 0
    for j in l:
        if j[0] in bset:
            return False
        #if 'R/B' in j[0]:
        #    return False
    return True


if 0:
    print "#!/bin/bash"
    print
    comm_strings = set()
    for info in info_triplets:
        req = from_time_get_space_track_request(info[3])

        if not os.path.exists(req[1]):
            comm_strings.add(req[0])
        
    for c in comm_strings:
        print """curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=NLHARR@berkeley.edu&password=asdfasdfasdfasdf'"""
        print c
    print """curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org//ajaxauth/logout"""

if 0:
    ov = []
    for info in info_triplets:
        print('**************************\n')
        print(info)
        
        sats = generate_satellite_list('3le_%s.tles'%t_to_date_str(info[3]))
        out_vals = get_observed_ra_dec(sats, info[3], info[1], info[2], do_parkes = False)

        ov.append(out_vals)
        print('\n**************************')

# generate event counts for various sig cutoffs
if 1:
    #['0 GRACE 2', '0 GRACE 1', '0 SWARM B', '0 SWARM A', '0 SWARM C'])
  import os.path, glob
  #ldf_fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_phinal_v2_G98_*pkl')
  ldf_fns = glob.glob('/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_phinal_v2_G98_ldfs_S3p0_2p0W4p5H0p8Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1.pkl')
  for ldf_fn in ldf_fns:


    #ldf_fn = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/damn_time_ldfs_S3p0W4p0H0p5Nv0p0Xv1p0E-04P4p0E-03Ng0Gr2Co6p0.pkl'
    
    #ldf_fn = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_ldfs_S3p0_2p0W4p0H0p5Nv0p0Xv1p5E-04P1p1E-03Ng0Gr2Co6p5Fs1.pkl'


    #ldf_fn = '/home/nlharr/spt3g_software/scratch/nlharr/frbhunt/ldfs/ldfs_phinal_ldfs_S3p0_2p0W4p0H0p5Nv0p0Xv1p5E-04P1p3E-03Ng0Gr2Co6p5Fs1.pkl'
    info_tag = os.path.basename(ldf_fn).split('.')[0].split('GLTCH')[-1][1:]

    #out_fn = 'sat_filtered_count_map_%s.pkl'%info_tag

    #out_fn = 'full_sat_filtered_thinner_%s.pkl'%info_tag
    out_fn = 'dummy.pkl'
    print(out_fn)

    if os.path.exists(out_fn):
        print('skipping', ldf_fn)
        continue


    sigs = range(7,17)
    count_map = {}

    import pickle
    evs = pickle.load(open(ldf_fn))[1][0].pixel_evs

    co_sig = 7.0
    evs = filter( lambda ev: min(ev.det_info[0].significance,ev.det_info[1].significance) > co_sig, evs)
    in_len = len(evs)

    pre_cm = {}
    for co_sig in sigs:
        pre_cm[co_sig] = len(filter( lambda ev: min(ev.det_info[0].significance,ev.det_info[1].significance) > co_sig, evs))

    out_evs = []
    out_names = []
    any_count = 0
    for i, ev in enumerate(evs):
        print(i, len(evs))
        ra = ev.det_info[0].ra
        ra = ra if ra > 0 else ra + 2 * np.pi
        
        info = (i, ra, ev.det_info[0].dec, ev.event_time.isoformat())
        sats = generate_satellite_list('3le_%s.tles'%t_to_date_str(info[3]))

        
        if False:
            out_vals = get_observed_ra_dec(sats, info[3], info[1], info[2], do_parkes = False)
        else:
            time = (core.G3Time(info[3]) - (np.random.randint(2) - 0.5) * 2 * (12 + np.random.rand()*12) *core.G3Units.hours)
            out_vals = get_observed_ra_dec(sats, time.isoformat(), info[1], info[2], do_parkes = False)    
            
        if len(out_vals) == 0:
            any_count += 1
        print("checking list")
        if check_lst(out_vals, do_all = False):
            out_evs.append(ev)
            ons = set()
            for ov in out_vals:
                ons.add(ov[0])
            out_names.append(list(ons))            
        print('^^^', ev.event_time.isoformat())
    print("Clean Count", any_count, len(evs))
    evs = out_evs
    for co_sig in sigs:
        count_map[co_sig] = len(filter( lambda ev: min(ev.det_info[0].significance,ev.det_info[1].significance) > co_sig, evs))

    print( count_map[sigs[0]] - len(filter(lambda x: x, out_names))  )

    print(count_map)
    print(pre_cm)
    
    
    if out_fn != 'dummy.pkl':
        import pickle
        pickle.dump(count_map, open(out_fn, 'w'))

'''
test_time = core.G3Time('2016-06-28T16:16:15.777592010')

print("generating sats")
sats = generate_satellite_list(test_time)

print("calculating ra decs")
ras, decs = get_observed_ra_dec(sats, test_time)
'''
