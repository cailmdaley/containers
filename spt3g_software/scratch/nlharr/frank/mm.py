from spt3g import core, mapmaker, coordinateutils, mapspectra
import spt3g.mapspectra.basicmaputils as BMU
import spt3g.mapspectra.modemixing as MM
import spt3g.mapspectra.fancyestimators as FE
import numpy as np
import pylab as pl
import pickle
from copy import copy



#apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/apod_mask.pkl'))
#hole_apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/apod_mask_hole.pkl'))

apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/wider_apod_mask.pkl'))
hole_apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/holey_wider_apod_mask.pkl'))
pnt_src_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/point_source_mask.pkl'))


ell_start = 70
ell_stop  = 90

cls = np.zeros((4,ell_stop))
cls[1, ell_start:ell_stop] = 1



#inpaint comparison
if 0:
    ell_bins = BMU.get_reg_spaced_ell_bins(range(10,3000,2))
    ells = map(np.mean, ell_bins)
    print('generating sky')
    t,q,u,cls = MM.generate_flat_sky(cls, apod_mask)
    q, u = BMU.flatten_pol(q,u)     

    print('estimating spectra')
    tt,ee,bb = MM.reference_power_spectra_function(ell_bins,t,q,u,apod_mask)

    no_in_chi_bb = FE.get_chi_b_spectra(ell_bins, q, u, hole_apod_mask)
    no_in_pure_bb = FE.get_pure_spectra(ell_bins, q, u, hole_apod_mask)

    in_chi_chi_bb = FE.get_chi_b_spectra(ell_bins, q, u, apod_mask, inpaint_mask=pnt_src_mask)
    
    print('hole spectra')
    htt,hee,hbb = MM.reference_power_spectra_function(ell_bins,t,q,u,hole_apod_mask)

    mapspectra.inpaint_map_laplace(pnt_src_mask, 10000, q)
    mapspectra.inpaint_map_laplace(pnt_src_mask, 10000, u)

    print('inpaint spectra')
    itt,iee,ibb = MM.reference_power_spectra_function(ell_bins,t,q,u,apod_mask)
    in_chi_bb = FE.get_chi_b_spectra(ell_bins, q, u, apod_mask)
    in_pure_bb = FE.get_pure_spectra(ell_bins, q, u, apod_mask)

    pf = pl.semilogy

    pf(ells, iee)
    pf(ells, hee)

    pf(ells, ibb)
    #pf(ells, hbb)

    pf(ells, in_chi_bb)
    pf(ells, in_pure_bb)
    
    pf(ells, no_in_chi_bb)
    pf(ells, no_in_pure_bb)
    pf(ells, in_chi_chi_bb)

    leg = pl.legend(['EE Inpaint', 'EE Holey',
                     'BB',  
                     '$\chi_B$ Inpaint', 
                     'Pure Inpaint', 
                     '$\chi_B$ Holey', 
                     'Pure Holey',
                     '$\chi_B$ Inpaint of $\chi_B$'
                 ])
    for line in leg.get_lines():
        line.set_linewidth(4.0)

    pl.xlabel('$l$')
    pl.ylabel('$C_{X,l}$')
    pl.title('Pseudo $C_l$ Esimator Comparison')
    pl.show()


#map projection comparison
if 1:

    #apod_mask = pickle.load(open('/home/nlharr/spt3g_software/scratch/nlharr/frank/tukey_apod_mask_0p25.pkl'))
    cs = ['r','g','b']
    ell_bins = BMU.get_reg_spaced_ell_bins(range(10,3000,10))
    ells = map(np.mean, ell_bins)

    desired_projections = [coordinateutils.MapProjection.ProjSansonFlamsteed,
                           coordinateutils.MapProjection.ProjStereographic,
                           coordinateutils.MapProjection.ProjLambertAzimuthalEqualArea]
    labels = ['Sanson Flamsteed', 
              'Stereographic', 
              'Lambert Azimuthal Equal Area']
    labels = map(lambda s: '$\chi_B$ '+s, labels)

    tts = []
    ees = []
    bbs = []

    chi_bbs = []
    for i, proj in enumerate(desired_projections):
        apod_mask_new = copy(apod_mask)
        apod_mask_new.proj = proj
        mapmaker.rebin_map(apod_mask, apod_mask_new)
        
        t,q,u,cls = MM.generate_flat_sky(cls, apod_mask_new)

        q, u = BMU.flatten_pol(q,u)
        tt,ee,bb = MM.reference_power_spectra_function(ell_bins,t,q,u,apod_mask_new)
        chi_bb = FE.get_chi_b_spectra(ell_bins, q, u, apod_mask)

        tts.append(tt)
        ees.append(ee)
        bbs.append(bb)
        chi_bbs.append(chi_bb)

    pl.plot(ells, ees[1])
    pl.plot(ells, chi_bbs[0])
    pl.plot(ells, chi_bbs[1])
    pl.semilogy(ells, chi_bbs[2])

    pl.legend(['EE Spectrum',labels[0],labels[1],labels[2]])
    pl.xlabel('$l$')
    pl.ylabel('$C_{X,l}$')
    pl.title('Flat Sky Projection Comparison')

    pl.show()


