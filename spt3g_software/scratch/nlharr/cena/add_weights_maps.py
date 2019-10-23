from spt3g import core, mapmaker, sptpol, todfilter, dfmux, util, timestreamflagging
import spt3g.mapmaker.mapmakerutils as mmu
import spt3g.mapmaker.summingmaps as SM
import spt3g.std_processing as std_processing
from copy import copy
import numpy as np
import pickle

cena_line_up = {
    'cena-20150312_131911.g3':'rcw38-20150312_123711.g3',
    'cena-20150325_151459.g3':'rcw38-20150325_143300.g3',
    'cena-20150301_230758.g3':'rcw38-20150301_222600.g3',
    'cena-20150314_025927.g3':'rcw38-20150314_021727.g3',
    'cena-20150303_053730.g3':'rcw38-20150303_045532.g3',
    'cena-20150315_161543.g3':'rcw38-20150315_153345.g3',
    'cena-20150303_122107.g3':'rcw38-20150303_113907.g3',
    'cena-20150317_065622.g3':'rcw38-20150317_061426.g3',
    'cena-20150304_174541.g3':'rcw38-20150304_170344.g3',
    'cena-20150318_194938.g3':'rcw38-20150318_190454.g3',
    'cena-20150401_014312.g3':'rcw38-20150401_103312.g3',
    'cena-20150306_071532.g3':'rcw38-20150306_063336.g3',
    'cena-20150320_091251.g3':'rcw38-20150320_083055.g3',
    'cena-20150307_204739.g3':'rcw38-20150307_200253.g3',
    'cena-20150321_225618.g3':'rcw38-20150321_221418.g3',
    'cena-20150309_094036.g3':'rcw38-20150309_085837.g3',
    'cena-20150322_114432.g3':'rcw38-20150322_110232.g3',
    'cena-20150310_235454.g3':'rcw38-20150310_231255.g3',
    'cena-20150324_013009.g3':'rcw38-20150324_004811.g3',
}

out_dir_2 = '/data/nlharr/cena_maps/pass_2/'
sum_dir = '/data/nlharr/cena_maps/map_sums/'

in_maps = map(lambda s:  out_dir_2 + 'map_' + s ,cena_line_up.keys())
point_source_fn = '/home/nlharr/spt_code/spt3g_software/scratch/nlharr/cena/cenamask.pkl'
mask = pickle.load(open(point_source_fn))

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = in_maps)
pipe.Add(SM.MakeMapsWeightedFromMap, mask = mask)
pipe.Add(SM.AddUnpolarizedMapsComingThrough)
pipe.Add(mapmaker.mapmakerutils.RemoveWeightModule)
pipe.Add(lambda fr: fr.type == core.G3FrameType.Map or 
         fr.type == core.G3FrameType.EndProcessing)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename=sum_dir+'summap.g3')
pipe.Run()
