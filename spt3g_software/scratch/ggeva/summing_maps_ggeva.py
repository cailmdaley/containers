# Modeled off of Daniel Dutcher's script
# spt3g_software/scratch/ddutcher/summing_maps.py

# How to run locally:
# python summing_maps_ggeva.py

from spt3g import core, coordinateutils, util
from glob import glob
import pickle as pkl

in_maps = sorted(glob('/spt/user/ggeva/PMNJ0210-5101-pixelraster/coadd_lr_0.25res_pol/*/coadd_lr_0.25res_*.g3'))
output_root = '/spt/user/ggeva/map_sums/PMNJ0210-5101-pixelraster/coadd_lr_0.25res_pol/'

# good=[]
# map_cuts = pkl.load(open('/home/ddutcher/data/bad_maps_5sig.pkl','rb'))
# for obs in in_maps:
#     obsid = (obs.split('/')[-1].split('_')[-1].replace('.g3',''))
#     if obsid not in map_cuts.keys():
#         good.append(obs)
#     else:
#         if 'med_wrms2' not in map_cuts[obsid]:
#             good.append(obs)

def foo(fr):
    if fr.type == core.G3FrameType.Observation:
        print(fr)

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader,filename = in_maps)
pipe.Add(foo)
# pipe.Add(util.framecombiner.MapFrameCombiner, fr_id='150GHz')
pipe.Add(util.framecombiner.MapFrameCombiner, fr_id='*Left-90GHz')
pipe.Add(util.framecombiner.MapFrameCombiner, fr_id='*Left-150GHz')
pipe.Add(util.framecombiner.MapFrameCombiner, fr_id='*Left-220GHz')
pipe.Add(util.framecombiner.MapFrameCombiner, fr_id='*Right-90GHz')
pipe.Add(util.framecombiner.MapFrameCombiner, fr_id='*Right-150GHz')
pipe.Add(util.framecombiner.MapFrameCombiner, fr_id='*Right-220GHz')
pipe.Add(lambda fr: 'Id' in fr and 'combined' in fr['Id'])
pipe.Add(core.G3Writer, filename=output_root + 'pmnj0210-5101-pixelraster_coadd_lr_0.25res_sum.g3')
pipe.Run()
