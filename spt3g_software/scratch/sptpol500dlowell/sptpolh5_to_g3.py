import glob, copy
import sptpol_software.util.hdf as hdf
import spt3g
from spt3g import core,dfmux,std_processing, mapmaker,  coordinateutils
import numpy as np
import argparse as ap

'''
This code takes sptpol map objects from sptpol_sofware, and put them into 
spt3g style map objects. It reads in sptpol .h5 map files and saves them as 
spt3g .g3 files. 

Example use : 

python spt3g_software/scratch/sptpol500dlowell/sptpolh5_to_g3.py /home/javva/spt3g_software/scratch/javva/bundles/bundles_090.h5 -o bundles_090_sptpol.g3

'''
P = ap.ArgumentParser(description='Sptpol maps into .g3 maps',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_file', action='store', nargs='+', default=[],
               help='Input files')
P.add_argument('-o', '--output', action='store', default='g3_from_pol_map.g3',
           help='Output 3g map filename')

args = P.parse_args()

if len(args.input_file) == 1:
    filename = args.input_file[0]
    data = hdf.readSptHDF5(filename)
    sm =  data.getStokesMatrixContents()# weighted map, shape (ny, nx, 1, [T,Q,U])
    smT = sm[::-1,::-1,0,0]
    smQ = sm[::-1,::-1,0,1]
    smU = sm[::-1,::-1,0,2]

    #instantiate map objects for each one, fill in meta data
    t_map = coordinateutils.FlatSkyMap(
        smT*core.G3Units.kelvin, 
        data.initArgs()['reso_arcmin']*core.G3Units.arcmin, 
        is_weighted = data.weighted_map, 
        proj = spt3g.coordinateutils.MapProjection(data.initArgs()['projection']),
        alpha_center = data.initArgs()['center'][0]*core.G3Units.degrees, 
        delta_center = data.initArgs()['center'][1]*core.G3Units.degrees, 
        coord_ref=spt3g.core.MapCoordReference.Equatorial, 
        units =  spt3g.core.G3TimestreamUnits.Tcmb, 
        pol_type=spt3g.core.MapPolType.T)

    q_map = copy.copy(t_map)
    np.asarray(q_map)[:,:] = smQ*core.G3Units.kelvin
    q_map.pol_type = spt3g.core.MapPolType.Q

    u_map = copy.copy(t_map)
    np.asarray(u_map)[:,:] = smU*core.G3Units.kelvin
    u_map.pol_type = spt3g.core.MapPolType.U

    #define a weights map
    WeightsMap = core.G3SkyMapWeights(t_map)
    sptpol_weights = data.initArgs()['weight']

    np.asarray(WeightsMap.TT)[:,:] = sptpol_weights[::-1,::-1,0,0]
    np.asarray(WeightsMap.TQ)[:,:] = sptpol_weights[::-1,::-1,1,0]
    np.asarray(WeightsMap.TU)[:,:] = sptpol_weights[::-1,::-1,2,0]
    np.asarray(WeightsMap.QQ)[:,:] = sptpol_weights[::-1,::-1,1,1]
    np.asarray(WeightsMap.QU)[:,:] = sptpol_weights[::-1,::-1,1,2]
    np.asarray(WeightsMap.UU)[:,:] = sptpol_weights[::-1,::-1,2,2]

    pipe = core.G3Pipeline()
    pipe.Add(core.G3InfiniteSource,type= spt3g.core.G3FrameType.Map,n=0)
    pipe.Add(mapmaker.MapInjector, map_id='sptpol_map',
             maps_lst=[t_map, q_map, u_map, WeightsMap], 
             is_stub=False, do_weight=True, make_polarized = True)
    pipe.Add(core.Dump)
    pipe.Add(core.G3Writer, filename = args.output)
    pipe.Run(profile=True)

else:
    print "Currently only supports changing one file at a time."
