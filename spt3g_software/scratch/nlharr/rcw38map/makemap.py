from spt3g import core, mapmaker
import spt3g.mapmaker.mapmakerutils as mmu
import spt3g.std_processing as std_processing


def doshit(dyno):
    fn = '/home/nlharr/tmp/rcw38-21Apr2015.g3'
    smstub = std_processing.CreateSourceMapStub('rcw38', x_len = 400, y_len = 400, 
                                                res = 0.5 * core.G3Units.arcmin,
                                                proj = mapmaker.MapProjection.ProjLambertAzimuthalEqualArea)
    individual_bolos_to_map = ['Sq1SBpol03Ch1', 'Sq1SBpol03Ch10', 
                               'Sq1SBpol03Ch11', 'Sq1SBpol03Ch12', 
                               'Sq1SBpol03Ch2', 'Sq1SBpol03Ch3', 
                               'Sq1SBpol03Ch4', 'Sq1SBpol03Ch5', 'Sq1SBpol03Ch6']
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader,
             filename = fn)
    pipe.Add(std_processing.MakeBootstrappedIndividualBoloMap,
             map_in = smstub, poly_order = 4,
             ts_in_key = 'RawTimestreams_I', 
             individual_bolos_to_map = individual_bolos_to_map,
             use_dyno_filter = dyno
    )
    map_extractor = mapmaker.mapmakerutils.ExtractTheMaps()
    pipe.Add(core.Dump)
    pipe.Add(map_extractor)
    pipe.Run(profile = True)
    return map_extractor

msked = doshit(True)
#reg = doshit(False)
