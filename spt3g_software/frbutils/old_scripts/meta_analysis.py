import numpy as np
import os.path
from glob import glob
import cPickle as pickle

def lazy_split_num(s):
    head = s.rstrip('0123456789')
    tail = s[len(head):]
    return head, tail

def parse_file_tag(fn):
    assert(os.path.basename(fn).split('_')[0] == 'info')
    assert(os.path.basename(fn).split('_')[1] == 'summary')
    #returns a dictionary mapping the key in teh label to its numerical value in a floating point number.
    #if broken, just rewrite
    return dict(map(lambda x: (x[0], float(x[1])), map(lazy_split_num, ''.join(os.path.basename(fn).split('.')[:-1]).split('_')[2:])))

def get_summary_stats(fn, has_injected):
    d = pickle.load(open(fn))
    info = {}
    #detection percentage of events
    for k in d:
        print k
        if has_injected:
            info[k] = {
                'injected_num_scans': d[k][0].n_scans_looked_through,
                'injected_num_found': d[k][0].n_events_found
                }
        else:
            info[k] = {}
        det_info = d[k][1]
        n_total_hits = 0
        for m in det_info.char_dic:
            for l in det_info.char_dic[m]:
                n_total_hits += sum(det_info.char_dic[m][l]['n_hits'])
        info[k]['n_total_hits'] = n_total_hits

    return info

info_files = glob('/data/spt/processed_data/*/info_summary_*.pkl')


out_lst = []
for info_file in info_files:
    fn_info = parse_file_tag(info_file)
    summary_stat = get_summary_stats(info_file, fn_info['injfs'])
    out_lst.append((fn_info, summary_stat))
