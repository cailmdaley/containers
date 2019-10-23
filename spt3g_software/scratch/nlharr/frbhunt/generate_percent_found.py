from spt3g import core, mapmaker, frbutils
import numpy as np
import pickle
from copy import copy

BY_SQUID = 1
BY_WAFER = 2
BY_NOT_SQUID = 3
BY_NOT_BOARD = 4

low_ll_cutoffs = [6,7,8, 9, 10, 11, 12, 13, 14, 15, 16]
#low_ll_cutoffs = [14]

def AddPixelIdSptpol(frame, bkey = 'BolometerProperties'):
    if frame.type == core.G3FrameType.Calibration:
        bprops = copy(frame[bkey])
        for k in bprops.keys():
            phys_name = bprops[k].physical_name
            bprops[k].pixel_id = phys_name[:phys_name.rfind('.')]
        del frame[bkey]
        frame[bkey] = bprops

@core.scan_func_cache_data(bolo_props = 'BolometerProperties')
def AddMapToZero(frame, store_key='ZeroMap', bolo_props = None):
    a = core.G3MapDouble()
    for k in bolo_props.keys():
        a[k] = 0
    frame[store_key] = a


class GetPf(object):
    def __init__(self):
        self.n_frames = 0
        self.ev_keys = ['FrbEvents',
                        'FrbEvsGl', 'FrbEvsGlCoSig', 'FrbEvsGlCoSigWaf', 'FrbEvsGlCoSigWafHev',
                        'FrbEvsGlCoSigWafHevNegEv', 'FrbEvsGlCoSigWafHevNegEvTwo', 
                        'FrbEvsGlCoSigWafHevNegEvTwoPnt', 'FrbEvsGlCoSigWafHevNegEvTwoPntSqs',
                        'FrbEvsGlCoSigWafHevNegEvTwoPntSqsWaf', 
                        #'FrbEvsGlCoSigWafHevNegEvTwoPntSqsWafIndS', 
                        #'FrbEvsGlCoSigWafHevNegEvTwoPntSqsWafIndSPha', 
                        'FrbEventsFilteringOut']
        self.n_found = {}
        for k in self.ev_keys:
            self.n_found[k] = 0
    def __call__(self, frame):
        if frame.type != core.G3FrameType.Scan:
            return

        #print(frame)


        self.n_frames += 1
        for k in self.ev_keys:
            evs = frame[k]
            for ev in evs:
                npair = 0
                for di in ev.det_info:
                    if di.was_injected:
                        npair += 1
                    if npair >= 2:
                        #print("Found", self.n_found, self.n_frames)
                        self.n_found[k] += 1



def generate_pf( pos_hist_bins, input_files, low_ll_cutoffs, 
                 grouping_type, output_prefix, 
                 input_key = 'FrbEvsGl',
                 squid_max_cutoff = 4, 
                 squid_mean_cutoff = 4, 
                 
                 wafer_ll_cutoff = 4, heavi_cutoff = 0.5,
                 min_var = 3e-5, max_var = 10e-5, percent_above_5_sig = 0.005,
                 geometry_index = 0, is_neg = False, bad_channels = None,
                 store_sig = -1, filt_evs_of_sign = False, co_sig = 7,
                 point_source_file = None, filt_same = True
                 ):

    output_file = (output_prefix + 
                    ('ldfs_S%.1f_%.1fW%.1fH%.1fNv%.1fXv%.1EP%.1ENg%dGr%dCo%.1fFs%d' % 
                     (squid_max_cutoff, squid_mean_cutoff, wafer_ll_cutoff, 
                      heavi_cutoff, min_var, max_var,
                      percent_above_5_sig, int(is_neg), grouping_type, 
                      co_sig, filt_same)).replace('.','p')
                   +'.pkl'
                   )
    print(output_file)

    pipe = core.G3Pipeline()
    pipe.Add( core.G3Reader, filename = input_files)#, fault_tolerant = True )
    pipe.Add(core.Delete, keys = ['FrbEvsGl'], type = core.G3FrameType.Scan)

    pipe.Add(AddMapToZero)
    pipe.Add(core.Dump, type = core.G3FrameType.Observation)

    pipe.Add(frbutils.frbfiltering.FrbFilterGlitchyDetectors,
             unlikely_cutoff = 0.99,
             input_frb_event_key = 'FrbEvents',
             output_frb_event_key = 'FrbEvsGl',
             n_tses_key = 'NumLiveDetectors')
    #pipe.Add(core.Dump)
    pipe.Add(frbutils.frbfiltering.SptpolFrbFilteringSecond,
             input_frb_event_key = input_key, 
             is_neg = is_neg,
             squid_max_cutoff = squid_max_cutoff,
             squid_mean_cutoff = squid_mean_cutoff,
             wafer_ll_cutoff = wafer_ll_cutoff,
             heavi_percent_of_amp = heavi_cutoff,
             filt_evs_of_sign = filt_evs_of_sign,
             point_source_file = point_source_file,
             filt_same = filt_same,
             co_sig=co_sig,
             del_big_ev_lists = False,
             phase_amp = None,
             )

    gpf = GetPf()
    pipe.Add(gpf)


    pipe.Add(frbutils.frbfiltering.ConstructValidBidsList,
             output_key = 'GoodChannels',
             var_min = min_var, var_max = max_var, 
             percent_above_5_sig = percent_above_5_sig,
             bad_channels = bad_channels)

    #pipe.Add(core.Dump, type = core.G3FrameType.Observation)
    ldfs = []
    for i, low_ll in enumerate(low_ll_cutoffs):
        sig_evs_name = 'FrbEvsSig%.1f' % low_ll
        pipe.Add(frbutils.frbfiltering.FilterSig,
                 event_key_in = 'FrbEventsFilteringOut',
                 event_key_out = sig_evs_name,
                 min_sig = low_ll,
                 max_sig = 1e12,
                 max_small_sig = None)


        force_not_shared_board = False
        force_not_shared_squid = False
        do_by_wafer = True
        do_by_squid = False
        if grouping_type == BY_SQUID:
            do_by_squid = True
            do_by_wafer = False
        elif grouping_type == BY_WAFER:
            pass
        elif grouping_type == BY_NOT_SQUID:
            force_not_shared_squid = True
        elif grouping_type == BY_NOT_BOARD:
            force_not_shared_board = True
        else:
            raise RuntimeError("Type no make sense")


        print("FNS ", force_not_shared_board, force_not_shared_squid)
        do_geometry = False
        store_sig = -1


        ldf = frbutils.frbanalysis.CosmicRayLdf(
            pos_hist_bins, sig_evs_name, 'GoodChannels',
            do_by_squid=do_by_squid, do_by_wafer=do_by_wafer,
            do_geometry=do_geometry,  
            get_pf = True,
            skip_pixel_evs = not (low_ll == store_sig),
            force_not_shared_board = force_not_shared_board,
            force_not_shared_squid = force_not_shared_squid)

        pipe.Add(ldf)
        ldfs.append(ldf)
    #pipe.Add(core.InjectDebug, type = core.G3FrameType.Scan)
    pipe.Run()

    pickle.dump( (low_ll_cutoffs, ldfs, gpf.n_frames, gpf.n_found), open(output_file, 'w') )




if __name__ == '__main__':
    from glob import glob
    from spt3g import mapmaker
    #input_files = sorted(glob('/spt/user/nlharr/output/search_66_wafer_tight_search/frb_search_66_wafer_ra23h30dec-55_idf_*_150ghz_processed.g3'))
    #input_files = sorted(glob('/spt/user/nlharr/output/search_66_wafer_2/*_processed.g3'))[::-1]
    point_source_file = '/spt/public/nlharr/frb_side_products/large_field_pntsrc_mask.fits'



    #for inj in ['4','8','11p6', '16','32','64','128','256'][::-1]:
    for inj in ['4','8','11p6'][::-1]:
    #for inj in ['64'][::-1]:
    #for inj in ['64'][::-1]:
        print("doing inj",inj)
        input_files = sorted(glob('/spt/user/nlharr/output/sq_inject_005_%s_66_wafer/*_processed.g3'%inj))
        #input_files = sorted(glob('/spt/user/nlharr/output/sq_inject_%s_66_wafer/*_processed.g3'%inj))

        #-39601430 -39755458
        #import pdb; pdb.set_trace()
        pos_hist_bins = np.arange(0,200,0.01)
        generate_pf( pos_hist_bins, input_files, low_ll_cutoffs, 
                     #grouping_type = BY_NOT_SQUID, 
                     #grouping_type = BY_NOT_BOARD, 
                     #grouping_type = BY_SQUID, 
                     grouping_type = BY_WAFER, 

                     output_prefix = 'pfs/pf_phinal_NEW_ALL_SQ_GLTCH_BE_G99_p75_%s_'%inj,
                     #output_prefix = 'pfs/pf_phinal_NEW_ALL_SQ_GLTCH_G99_p80_%s_'%inj,
                     #output_prefix = 'pfs/pf_phinal_005Width_NEW_ALL_SQ_GLTCH_G99_p80_%s_'%inj,
                     #output_prefix = 'pfs/phase_test_p1_%s_'%inj,


                     squid_max_cutoff = 3.0, 
                     squid_mean_cutoff = 2.0,

                     wafer_ll_cutoff = 5.0, 
                     min_var = 2.0e-5,
                     max_var = 15e-5,
                     
                     heavi_cutoff = 0.75,
                     
                     percent_above_5_sig = 0.0011, 

                     bad_channels = ['Sq6SBpol05Ch5', 'Sq6SBpol05Ch9'],                     
                     #is_neg = False,
                     is_neg = True,
                     store_sig = 6,
                     filt_evs_of_sign = True, co_sig = 6.5,
                     point_source_file  = point_source_file,
                     filt_same = True
                     )
