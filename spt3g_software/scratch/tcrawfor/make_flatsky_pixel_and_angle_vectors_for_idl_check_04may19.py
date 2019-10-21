import numpy as np
import scipy
import pickle
from spt3g import core, std_processing, gcp, coordinateutils, mapmaker

# using 

#outdir = '/spt/user/tcrawfor/public/pix_ang_tests/'
outdir = '/spt/user/tcrawfor/public/pix_ang_tests_2/'

x_len = 512
y_len = 256
nrand = 10000
reso_arcmin = 1.0
res = reso_arcmin*core.G3Units.arcmin
projections = ['ProjSansonFlamsteed','ProjLambertAzimuthalEqualArea']#,'ProjBICEP']
ra = 12.783*core.G3Units.deg
dec = -41.871*core.G3Units.deg # sort of random numbers
ddec = y_len*res
dra = x_len*res*np.cos(dec-ddec/2.)


for projection in projections:
    proj = coordinateutils.MapProjection.names[projection]
    print(proj)
    map_params = coordinateutils.FlatSkyMap(
        x_len = x_len, y_len = y_len, res = res,
        proj = proj, alpha_center = ra, delta_center = dec,
        pol_type = core.MapPolType.T,
        coord_ref = core.MapCoordReference.Equatorial)
    iptemp_in = np.arange(y_len*x_len)
    iptemp_in.tofile(outdir + 'iptemp_in_' + projection)
    radec_out = np.asarray([map_params.pixel_to_angle(ip) for ip in iptemp_in])
    radec_out.tofile(outdir + 'radec_out_' + projection)
    xtemp_in = np.mod(iptemp_in,x_len)
    ytemp_in = iptemp_in/x_len
    xtemp_in_2 = xtemp_in + np.random.rand(len(xtemp_in)) - 0.5
    ytemp_in_2 = ytemp_in + np.random.rand(len(ytemp_in)) - 0.5
    xtemp_in_2.tofile(outdir + 'xtemp_in_2_' + projection)
    ytemp_in_2.tofile(outdir + 'ytemp_in_2_' + projection)
    radec_out_2 = np.asarray([map_params.xy_to_angle(xt,yt) for xt,yt in zip(xtemp_in_2,ytemp_in_2)])
    radec_out_2.tofile(outdir + 'radec_out_2_' + projection)
    ra_rand = (np.random.rand(nrand)*dra - dra/2.)*0.9 + ra
    dec_rand = (np.random.rand(nrand)*ddec - ddec/2.)*0.9 + dec
    ra_rand.tofile(outdir + 'ra_rand_' + projection)
    dec_rand.tofile(outdir + 'dec_rand_' + projection)
    iptemp_out = np.asarray([map_params.angle_to_pixel(rt,dt) for rt,dt in zip(ra_rand,dec_rand)])
    iptemp_out.tofile(outdir + 'iptemp_out_' + projection)







