import glob
import sptpol_software.util.hdf as hdf
import spt3g
from spt3g import core,dfmux,std_processing, mapmaker,  coordinateutils
import numpy as np
#for filename in glob.glob('dictionary/*'):
for i in [1]:
    filename = 'spt3g/spt3g_software/scratch/javva/bundles/bundles_090.h5'
    data = hdf.readSptHDF5(filename)
    sm =  data.getStokesMatrixContents()# weighted map, shape (ny, nx, 1, [T,Q,U])
    smT = sm[:,:,0,0]
    smQ = sm[:,:,0,1]
    smU = sm[:,:,0,2]

    #instantiate map objects for each one, fill in meta data
    t_map = coordinateutils.FlatSkyMap(smT*core.G3Units.kelvin, data.initArgs()['reso_arcmin']*core.G3Units.arcmin, is_weighted = data.weighted_map, proj = spt3g.coordinateutils.MapProjection(data.initArgs()['projection']),alpha_center = data.initArgs()['center'][0], delta_center = data.initArgs()['center'][1], coord_ref=spt3g.core.MapCoordReference.Equatorial, units =  spt3g.core.G3TimestreamUnits.Kcmb, pol_type=spt3g.core.MapPolType.T)
    q_map = coordinateutils.FlatSkyMap(smQ*core.G3Units.kelvin, data.initArgs()['reso_arcmin']*core.G3Units.arcmin, is_weighted = data.weighted_map, proj = spt3g.coordinateutils.MapProjection(data.initArgs()['projection']),alpha_center = data.initArgs()['center'][0], delta_center = data.initArgs()['center'][1], coord_ref=spt3g.core.MapCoordReference.Equatorial, units =  spt3g.core.G3TimestreamUnits.Kcmb, pol_type=spt3g.core.MapPolType.Q)
    u_map = coordinateutils.FlatSkyMap(smU*core.G3Units.kelvin, data.initArgs()['reso_arcmin']*core.G3Units.arcmin, is_weighted = data.weighted_map, proj = spt3g.coordinateutils.MapProjection(data.initArgs()['projection']),alpha_center = data.initArgs()['center'][0], delta_center = data.initArgs()['center'][1], coord_ref=spt3g.core.MapCoordReference.Equatorial, units =  spt3g.core.G3TimestreamUnits.Kcmb, pol_type=spt3g.core.MapPolType.U)
    
    #define a weights map

    WeightsMap = core.G3SkyMapWeights(t_map)
    sptpol_weights = data.initArgs()['weight']

    #Fill in WeightsMap TT
    tt_w = sptpol_weights[:,:,0,0].copy(order='C') #makes contiguous so b.p. works with python slices
    tt__w = coordinateutils.FlatSkyMap(tt_w,data.initArgs()['reso_arcmin']*core.G3Units.arcmin)
    WeightsMap.TT = tt__w

    #Fill in WeightsMap TQ
    tq_w = sptpol_weights[:,:,0,1].copy(order='C') #makes contiguous so b.p. works with python slices
    tq__w = coordinateutils.FlatSkyMap(tq_w,data.initArgs()['reso_arcmin']*core.G3Units.arcmin)
    WeightsMap.TQ = tq__w

    #Fill in WeightsMap TU
    tu_w = sptpol_weights[:,:,0,2].copy(order='C') #makes contiguous so b.p. works with python slices
    tu__w = coordinateutils.FlatSkyMap(tu_w,data.initArgs()['reso_arcmin']*core.G3Units.arcmin)
    WeightsMap.TU = tu__w

    #Fill in WeightsMap QQ
    qq_w = sptpol_weights[:,:,1,1].copy(order='C') #makes contiguous so b.p. works with python slices
    qq__w = coordinateutils.FlatSkyMap(qq_w,data.initArgs()['reso_arcmin']*core.G3Units.arcmin)
    WeightsMap.QQ = qq__w

    #Fill in WeightsMap QU
    qu_w = sptpol_weights[:,:,1,2].copy(order='C') #makes contiguous so b.p. works with python slices
    qu__w = coordinateutils.FlatSkyMap(qu_w,data.initArgs()['reso_arcmin']*core.G3Units.arcmin)
    WeightsMap.QU = qu__w

    #Fill in WeightsMap UU
    uu_w = sptpol_weights[:,:,2,2].copy(order='C') #makes contiguous so b.p. works with python slices
    uu__w = coordinateutils.FlatSkyMap(uu_w,data.initArgs()['reso_arcmin']*core.G3Units.arcmin)
    WeightsMap.UU = uu__w

    pipe = core.G3Pipeline()
    pipe.Add(core.G3InfiniteSource,type= spt3g.core.G3FrameType.Map,n=0)
    pipe.Add(mapmaker.MapInjector, map_id='sptpol_map',maps_lst=[t_map, q_map, u_map, WeightsMap], is_stub=False, do_weight=True, make_polarized = True)
    pipe.Add(core.G3Writer, filename = 'bundles_%s_jt.g3'%k)
    pipe.Run(profile=True)

