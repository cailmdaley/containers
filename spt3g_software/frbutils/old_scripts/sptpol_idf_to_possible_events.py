import argparse
from spt3g import frbutils 

parser = argparse.ArgumentParser()

parser.add_argument('--input_file', default = '/home/nlharr/tmp/ra23h30dec-55_idf_20120707_163808_090ghz.h5')
parser.add_argument('--output_file', default = 'dummy.g3')
parser.add_argument('--model_width', type=int, default = 13)
parser.add_argument('--find_ll_thresh', type=float, default = 10.0)
parser.add_argument('--other_ll_thresh', type=float, default = 8.0)

parser.add_argument('--filter_poly_order', type=int, default = 11)
parser.add_argument('--fit_poly_order', type=int, default = 1)
parser.add_argument('--min_peak_distance', type=int, default = 30)

parser.add_argument('--inject_fake_signal', type=int, default = 0)
parser.add_argument('--injection_signficance', type=float, default = -1)

args = parser.parse_args()

if args.inject_fake_signal:
    assert(args.injection_signficance > 0)

frbutils.frbhunting.do_frb_search(idf_fn = args.input_file, 
                                  out_fn = args.output_file,
                                  
                                  filter_poly_order = args.filter_poly_order,
                                  fit_poly_order = args.fit_poly_order,
                                  model_width = args.model_width,
                                  find_ll_thresh = args.find_ll_thresh,
                                  other_ll_thresh = args.other_ll_thresh,
                                  min_peak_distance = args.min_peak_distance,
                                  
                                  inject_fake_signal = args.inject_fake_signal, 
                                  inject_significance = args.injection_signficance)
