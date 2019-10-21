from spt3g import core
from spt3g.frbutils import frbanalysis as frbanal


frame = frbanal.make_fake_scan_frame_from_frb_ev_lst(1, 'empty', 
                                                     hwm_folder = '/home/nlharr/tmp/frb_side_products/hwm_20140318/',
                                                     pkl_file = '/home/nlharr/tmp/frb_side_products/TES_positions_150ghz.pkl')


pairs = frbanal.try_to_find_close_pairs(frame)

print pairs
print map(lambda p: frbanal.get_bolo_pair_dist(frame, p[0], p[1]), pairs)
