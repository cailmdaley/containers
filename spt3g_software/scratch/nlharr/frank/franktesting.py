from spt3g import core, mapmaker, coordinateutils, mapspectra, mapmaker
import spt3g.mapspectra.basicmaputils as BMU
import spt3g.mapspectra.modemixing as MM
import spt3g.mapspectra.fancyestimators as FE
import numpy as np
import pylab as pl
import pickle
from copy import copy


def unweight_map_fake(t, q, u, pol_ang):
    sc = mapmaker.Mat1x3()
    mapmaker.set_stokes_coupling(pol_ang, 1, 1, sc)
    t = t.Clone(True)
    q = q.Clone(True)
    u = u.Clone(True)

    return t * sc[0] + q * sc[1] + u * sc[2]


def add_weight(s, pol_ang):
    sc = mapmaker.Mat1x3()
    mapmaker.set_stokes_coupling(pol_ang, 1, 1, sc)
    t = sc[0]* s
    q = sc[1]* s
    u = sc[2]* s
    w = core.G3SkyMapWeights(t)
    np.asarray(w.TT)[:] = sc[0]*sc[0]
    np.asarray(w.TQ)[:] = sc[0]*sc[1]
    np.asarray(w.TU)[:] = sc[0]*sc[2]
    np.asarray(w.QQ)[:] = sc[1]*sc[1]
    np.asarray(w.QU)[:] = sc[1]*sc[2]
    np.asarray(w.UU)[:] = sc[2]*sc[2]
    return t, q, u, w

def poly_filter_map(m, deg = 15):
    m = copy(m)
    nrows, ncols = np.shape(m)
    x = np.arange(ncols)
    for row in range(nrows):
        p = np.polyfit(x, np.asarray(m)[row,:], deg)
        np.asarray(m)[row,:] -= np.polyval(p, x)
    return m


def poly_filter_map_non_zero(m, deg = 15):
    m = copy(m)
    nrows, ncols = np.shape(m)
    x = np.arange(ncols)
    for row in range(nrows):
        y = np.asarray(m)[row,:]
        inds = np.where(y != 0)
        xs = x[inds]
        ys = y[inds]
        if len(ys) == 0:
            continue
        p = np.polyfit(xs, ys, deg)
        ys -= np.polyval(p, xs)
        y[inds] = ys
        np.asarray(m)[row,:] = y
    return m

def do_thing(t,q,u, pol_ang,p=5):
    s = unweight_map_fake(t, q, u, pol_ang)
    s = poly_filter_map_non_zero(s,p)
    return add_weight(s, pol_ang)

def do_other_thing(t,q,u, pol_ang,p=15):
    pol_angs = np.pi * np.array([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75])
    tsum,qsum,usum,wsum = do_thing(t, q, u, pol_angs[0], p)
    for pol_ang in pol_angs[1:]:
        ttmp,qtmp,utmp,wtmp = do_thing(t, q, u, pol_ang, p)
        tsum += ttmp
        qsum += qtmp
        usum += utmp
        wsum += wtmp
    tr, qr, ur = mapmaker.remove_weight(tsum, qsum, usum, wsum)
    return tr,qr,ur

def read_camb_cls(fn):
    f = open(fn)
    ells = [0, 1]
    tt = [0, 0]
    ee = [0, 0]
    bb = [0, 0]
    te = [0, 0]
    for line in f:
        if line.strip()[0] == '#':
            continue
        ls = line.strip().split()
        ells.append(int(ls[0]))
        tt.append(float(ls[1]))
        ee.append(float(ls[2]))
        bb.append(float(ls[3]))
        te.append(float(ls[4]))
    out_arr = np.zeros( (4, len(tt)))
    out_arr[0,:] = np.array(tt)
    out_arr[1,:] = np.array(ee)
    out_arr[2,:] = np.array(bb)
    out_arr[3,:] = np.array(te)
    return ells, out_arr



#apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/wider_apod_mask.pkl')) #good
#apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/tukey_apod_mask_0p25.pkl')) #orig


apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/frank_tukey.pkl')) 


apod_mask_new = copy(apod_mask)
apod_mask_new.proj = coordinateutils.MapProjection.ProjStereographic
mapmaker.rebin_map(apod_mask, apod_mask_new)
apod_mask = apod_mask_new

#inner_apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/widest_inner_apod_mask.pkl'))
#inner_apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/ridic_apod_mask.pkl'))
#inner_apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/inner_ridic_apod_mask.pkl')) #good
#inner_apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/tukey25_inner.pkl'))# orig

inner_apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/frank_tukey_inner_less.pkl'))
inner_apod_mask_new = copy(inner_apod_mask)
inner_apod_mask_new.proj = coordinateutils.MapProjection.ProjStereographic
mapmaker.rebin_map(inner_apod_mask, inner_apod_mask_new)
inner_apod_mask = inner_apod_mask_new

print('reading')
ells, cls = read_camb_cls('test_lensedCls.dat')

if 0:
    print('making maps')
    t,q,u,out_cls = MM.generate_flat_sky(cls, apod_mask)
    q, u = BMU.flatten_pol(q,u)     
    pickle.dump( (t,q,u), open('temp_frank_map.pkl', 'w'))
else:
    t,q,u = pickle.load(open('temp_frank_map.pkl'))

ell_bins = BMU.get_reg_spaced_ell_bins(range(10,3000,10))
x = map(np.mean, ell_bins)
fancy_qu_to_e_only_qu = FE.fancy_qu_to_e_only_qu

tr, qn, un = do_other_thing(t,q,u, 1)
qdebf, udebf, d = fancy_qu_to_e_only_qu(qn, un, apod_mask)

bad_inds = np.where(np.asarray(qdebf) == 0) 
np.asarray(t)[bad_inds] = 0
np.asarray(q)[bad_inds] = 0
np.asarray(u)[bad_inds] = 0

tdeb, qdeb, udeb = do_other_thing(t,qdebf,udebf, 4)
ts, qs, us = do_other_thing(t,q,u, 4)
qr = qs - qdeb
ur = us - udeb


tt,ee,bb = MM.reference_power_spectra_function(ell_bins,ts,qs,us,inner_apod_mask)

ftt,fee,fbb = MM.reference_power_spectra_function(ell_bins,t,qr,ur,inner_apod_mask)
fno_in_chi_bb = FE.get_chi_b_spectra(ell_bins, qr, ur, inner_apod_mask)
full_no_in_chi_bb = FE.get_chi_b_spectra(ell_bins, qs, us, inner_apod_mask)


pl.plot(x,tt)
pl.loglog(ells, cls[0,:])
pl.plot(x,ee)
pl.loglog(ells, cls[1,:])

pl.plot(x,fbb)
pl.plot(x,fno_in_chi_bb)
pl.plot(x,full_no_in_chi_bb)
pl.loglog(ells, cls[2,:])
pl.ion()

pl.legend(['TT', 'TT Theory', 'EE', 'EE Theory', 
           'FRANK $\chi_B$',  'FRANK BB', 'Poly $\chi_B$', 'Theory BB'])

pl.xlabel('$l$')
pl.ylabel('$C_{X,l}$')
pl.title('BB Estimation with FRANK, Toy Example')
pl.show()


if 0:
    tt,ee,bb = MM.reference_power_spectra_function(ell_bins,t,q,u,apod_mask)
    no_in_chi_bb = FE.get_chi_b_spectra(ell_bins, q, u, apod_mask)
    no_in_pure_bb = FE.get_pure_spectra(ell_bins, q, u, apod_mask)

if 0:
    fno_in_chi_bb = FE.get_chi_b_spectra(ell_bins, qr, ur, inner_apod_mask)
    fno_in_pure_bb = FE.get_pure_spectra(ell_bins, qr, ur, inner_apod_mask)

#pl.plot(x, tt)
#pl.plot(x,ee)
#pl.loglog(x,bb)
#pl.loglog(x,no_in_chi_bb)
#pl.loglog(x,ftt)
#pl.loglog(x,fno_in_chi_bb)
#pl.loglog(ells, cls[2,:])
#pl.show()
