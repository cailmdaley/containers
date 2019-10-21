import argparse
from spt3g import core, frbutils 


filter_pipelines = {'first_pass': (frbutils.frbfiltering.firstpass_filtering, True, True)
                    #'group_filter': (frbutils.frbfiltering.amp_group_filter, False, False)
                    }



parser = argparse.ArgumentParser()

parser.add_argument('--input_file')
parser.add_argument('--output_file')
parser.add_argument('--tes_pos_pkl')
parser.add_argument('--hwm_dir')
parser.add_argument('--filter_pipeline')

args = parser.parse_args()


assert(args.filter_pipeline in filter_pipelines)
finfo = filter_pipelines[args.filter_pipeline]
kwargs = {}
if finfo[1]:
    kwargs['hwm_fn'] = args.hwm_dir
if finfo[2]:
    kwargs['tes_pos_pkl_fn'] = args.tes_pos_pkl

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = args.input_file)
pipe.Add(finfo[0], **kwargs)
pipe.Add(core.G3Writer, filename = args.output_file)
pipe.Run()
