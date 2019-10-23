import matplotlib as mpl
mpl.use('Agg')
import numpy as np
from matplotlib import pylab as plt
pl=plt
from spt3g import core, std_processing
import sys
sys.path.append('/home/axf295/2019/code/spt3g_software/polcal/python/')
import CenA_Map_Utils as CMU
import General_Utils  as GU
import Circle_Functions as CF


BANDS    = ['90','150','220']
BANDFIG    = {'90':1,'150':2,'220':3}
BANDCOLORS = {'90': 'b','150': 'g','220':'r'}

def plot_pernomang_coadds(nomang_coadds,obsID,source,savedir='./'):
    bandfig = BANDFIG
    for wafer in nomang_coadds:
        numangs = len(nomang_coadds[wafer].keys())
        plt.figure(figsize=(15,5*numangs))
        angnum = 1
        for ang in sorted(nomang_coadds[wafer].keys()):
            for band in nomang_coadds[wafer][ang]:
                plt.subplot(numangs,3,3*(angnum-1)+bandfig[band])
                imdat = nomang_coadds[wafer][ang][band]
                plt.imshow(imdat)
                plt.colorbar(fraction=.046,pad=.04)
                plt.title('%sGHz, %s %s deg'%(band,wafer,ang))
            angnum+=1
        plt.suptitle('%s %s Nomang Coadds Obs %s'%(source,wafer,obsID))
        plt.savefig(savedir+'%s_%s_NomAng_Coaddmaps_%s.png'%(source,wafer,obsID))
        plt.close('all')
    return



def plot_perband_coadds(perbandcoadds,numperband={},obs = '',source='CenA',savedir='./'):
    bandfig=BANDFIG
    plt.figure(figsize=(15,5))
    for band in perbandcoadds:
            plt.subplot(1,3,bandfig[band])
            imdat = perbandcoadds[band]
            plt.imshow(imdat)
            plt.colorbar(fraction=.046,pad=.04)
            plt.grid(color='w',ls='--')
            if band in numperband:
                plt.title('%sGHz T Coadd: %s bolos'%(band,int(numperband[band])))
            else:
                plt.title('%sGHz T Coadd'%band)
    plt.suptitle('%s Per-Band Coadds %s'%(source,obs))
    plt.savefig(savedir+'%s_PerBand_Coaddmaps_%s.png'%(source,obs))
    plt.close('all')
    return


def hist_singlebolo_polangle_residuals(singlebolo_params,boloprops,savedir='./'):
    
    residuals      = {}
    res_sigma      = {}
    bestfit_angles = {}
    unique_nomangs = []
    for b in singlebolo_params:
        if b not in boloprops:
            continue
        
        band = str(int(boloprops[b].band/core.G3Units.GHz))
        if band == '-1':
            continue
        
        wafer = boloprops[b].physical_name.split('_')[0].upper()
        nomang = np.floor(np.rad2deg(boloprops[b].pol_angle))
        if nomang<0:
            nomang = 180+nomang
        elif nomang > 180:
            nomang = nomang-180
        
        if nomang not in unique_nomangs:
            unique_nomangs.append(nomang)
            
        if wafer not in residuals:
            residuals[wafer]      = {}
            res_sigma[wafer]      = {}
            bestfit_angles[wafer] = {}
        if nomang not in residuals[wafer]:
            residuals[wafer][nomang]      = {}
            res_sigma[wafer][nomang]      = {}
            bestfit_angles[wafer][nomang] = {}
        if band not in residuals[wafer][nomang]:
            residuals[wafer][nomang][band]      = []
            res_sigma[wafer][nomang][band]      = []
            bestfit_angles[wafer][nomang][band] = []
        
        residuals[wafer][nomang][band].append(singlebolo_params[b]['res'])
        res_sigma[wafer][nomang][band].append(singlebolo_params[b]['rawPhi'][1])
        phi = singlebolo_params[b]['rotPhi']
        
        bestfit_angles[wafer][nomang][band].append(phi)
    
    nomang_residuals = {}
    nomang_fits      = {}
    nomang_sigmas    = {}
    for ang in sorted(unique_nomangs):
        nomang_residuals[ang] = {}
        nomang_fits[ang]      = {}
        nomang_sigmas[ang]    = {}
        for band in BANDCOLORS:
            nomang_residuals[ang][band] = []
            nomang_fits[ang][band]      = []
            nomang_sigmas[ang][band]    = []
    ## Per-Wafer Per-band Per-Nominal Angle Plots
    for wafer in residuals:
        pl.figure(figsize=(20,10))
        angnum = 1
        for nomang in residuals[wafer]:
            pl.subplot(2,4,angnum)
            for band in residuals[wafer][nomang]:
                nomang_residuals[nomang][band] +=  residuals[wafer][nomang][band]
                nomang_sigmas[nomang][band] += res_sigma[wafer][nomang][band]
                pl.hist( residuals[wafer][nomang][band] ,bins=np.linspace(-90,90,61),
                         histtype='step', color=BANDCOLORS[band],label=band+'GHz')
            pl.grid()
            pl.title('%s Deg'%nomang)
            pl.xlabel('Residual Angle (deg)')
            pl.legend(loc='best')
            
            pl.subplot(2,4,angnum+4)
            for band in residuals[wafer][nomang]:
                nomang_fits[nomang][band] += bestfit_angles[wafer][nomang][band]
                pl.hist( bestfit_angles[wafer][nomang][band] ,bins=np.linspace(0,180,61),
                         histtype='step', color=BANDCOLORS[band],label=band+'GHz')
            pl.grid()
            pl.title('%s Deg'%nomang)
            pl.xlabel('Best Fit Angle (deg)')
            pl.legend(loc='best')
            
            angnum+=1
            
        pl.suptitle('%s Pol Angle Residuals and Best-Fits'%wafer)
        pl.savefig(savedir+'%s_PolAngle_Residuals_and_BestFits.png'%wafer)
        pl.close()
        
        
    perbandresiduals = {}
    perbandsigmas    = {}
    pl.figure(figsize=(20,20))
    angnum = 0
    for nomang in sorted(nomang_residuals.keys()):
        angnum+=1
        pl.subplot(4,4,angnum)
        pl.title('%s deg'%nomang)
        for band in nomang_residuals[nomang]:
            if band not in perbandresiduals:
                perbandresiduals[band] = []
                perbandsigmas[band]    = []
            perbandresiduals[band] += nomang_residuals[nomang][band]
            perbandsigmas[band]    += nomang_sigmas[nomang][band]
            m,s = GU.calc_median_and_1sigma(nomang_residuals[nomang][band],sigmaclip=3)
            pl.hist(GU.remove_outliers(nomang_residuals[nomang][band]),bins=np.linspace(-90,90,61),
                    histtype='step',color=BANDCOLORS[band],label='%s GHz : %.2f +- %.2f'%(band,m,s))
        pl.grid()
        pl.legend()
        pl.xlabel('Residual Angle (deg)')
    pl.suptitle('Per-Nomang Pol Angle Residuals and Best-Fits')
    pl.savefig(savedir+'PerNomang_PolAngle_Residuals_and_BestFits.png')
    pl.close()  
    
    pl.figure(figsize=(10,10))
    for band in perbandresiduals:
        m,s = GU.calc_median_and_1sigma(perbandresiduals[band],sigmaclip=3)
        pl.hist(perbandresiduals[band],bins=np.linspace(-90,90,61),
                histtype='step',color=BANDCOLORS[band],label='%s GHz : %.2f +- %.2f'%(band,m,s))
    pl.suptitle('Per-Band Residuals')
    pl.legend()
    pl.grid()
    pl.savefig(savedir+'PerBand_PolAngle_Residuals.png')
    pl.close('all')
     
    pl.figure(figsize=(10,10))
    for band in perbandsigmas:
        m,s = GU.calc_median_and_1sigma(perbandsigmas[band],sigmaclip=3)
        print(band,m,s)
        goodvals = np.copy(perbandsigmas[band])[np.where(np.isfinite(perbandsigmas[band]))]
        
        pl.hist(goodvals,bins=np.linspace(0,90,31),
                histtype='step',color=BANDCOLORS[band],label='%s GHz : %.2f +- %.2f'%(band,m,s))
    pl.suptitle('Per-Band 1-sigma Uncert in Fit')
    pl.legend()
    pl.grid()
    pl.savefig(savedir+'PerBand_PolAngle_Residual_Uncert.png')
    pl.close('all')
        
        
    return
            
            
def plot_coadd_fit_params(coaddparams,perwafer=False,savedir = './'):
    bandfig=BANDFIG
    ## Plot rotated best fit angles
    if not perwafer:
        pl.figure(figsize=(10,10))
    
    for wafer in coaddparams:
        if perwafer:
            pl.figure(figsize=(10,10))
        if wafer == 'all' :
            m = 's'
            a = .75
        elif not perwafer:
            a = .25
            m = 'o'
        else:
            a = .75
            m = 'o'
        for ang in coaddparams[wafer]:
            for band in coaddparams[wafer][ang]:
                pl.subplot(2,1,1)
                pl.errorbar(coaddparams[wafer][ang][band]['nom'],
                            coaddparams[wafer][ang][band]['rotPhi'],
                            yerr=coaddparams[wafer][ang][band]['rawPhi'][1],
                            marker=m,color=BANDCOLORS[band],alpha=a)
                
                pl.subplot(2,1,2)
                pl.errorbar(coaddparams[wafer][ang][band]['nom'],
                            coaddparams[wafer][ang][band]['res'],
                            yerr=coaddparams[wafer][ang][band]['rawPhi'][1],
                            marker=m,color=BANDCOLORS[band],alpha=a)
                
        if perwafer:
            pl.subplot(211)
            pl.plot([0,180],[0,180],color='k',ls='--')
            pl.xlabel('Nominal Angle [deg]')
            pl.ylabel('Best Fit angle [deg]')
            pl.ylim([-20,200])
            pl.xlim([-5,180])
            pl.grid()

            pl.subplot(212)
            pl.plot([0,180],[0,0],color='k',ls='--')
            pl.xlabel('Nominal Angle [deg]')
            pl.ylabel('Fit Residual [deg]')
            pl.ylim([-25,25])
            pl.grid()
            
            pl.suptitle('%s per band NomAngCoadd fits'%wafer)
            pl.savefig(savedir+'%s_perband_NomAngCoadd_fits.png'%wafer)
            pl.close('all')
    if not perwafer:
        pl.subplot(211)
        pl.plot([0,180],[0,180],color='k',ls='--')
        pl.xlabel('Nominal Angle [deg]')
        pl.ylabel('Best Fit angle [deg]')
        pl.ylim([-20,200])
        pl.xlim([-5,180])
        pl.grid()

        pl.subplot(212)
        pl.plot([0,180],[0,0],color='k',ls='--')
        pl.xlabel('Nominal Angle [deg]')
        pl.ylabel('Fit Residual [deg]')
        pl.ylim([-25,25])
        pl.grid()
        pl.suptitle('Per-Wafer per band NomAngCoadd fits')
        pl.savefig(savedir+'Perwafer_perband_NomAngCoadd_fits.png')
        pl.close('all')
    
    ## Plot raw fit angles
    pl.figure(figsize=(10,10))
    for wafer in coaddparams:
        if wafer == 'all':
            m = 's'
            a = 1.
        else:
            a = .25
            m = 'o'
        for ang in coaddparams[wafer]:
            for band in coaddparams[wafer][ang]:
                pl.subplot(2,1,1)
                pl.errorbar(coaddparams[wafer][ang][band]['nom'],
                            coaddparams[wafer][ang][band]['rawPhi'][0],
                            yerr=coaddparams[wafer][ang][band]['rawPhi'][1],
                            marker=m,color=BANDCOLORS[band],alpha=a)
                
                pl.subplot(2,1,2)
                pl.errorbar(coaddparams[wafer][ang][band]['nom'],
                            np.asarray(coaddparams[wafer][ang][band]['rawPhi'][0])-coaddparams[wafer][ang][band]['nom'],
                            yerr=coaddparams[wafer][ang][band]['rawPhi'][1],
                            marker=m,color=BANDCOLORS[band],alpha=a)
    pl.subplot(211)
    pl.plot([0,180],[0,180],color='k',ls='--')
    pl.xlabel('Nominal Angle [deg]')
    pl.ylabel('Raw Fit angle [deg]')
    pl.grid()
    
    pl.subplot(212)
    pl.plot([0,180],[0,0],color='k',ls='--')
    pl.xlabel('Nominal Angle [deg]')
    pl.ylabel('Raw Fit Residual [deg]')
    pl.grid()
    pl.suptitle('Per-Wafer per-band NomAngCoadd raw fits')
    pl.savefig(savedir+'Perwafer_perband_NomAngCoadd_rawfits.png')
    pl.close('all')
    
    
    return





def hist_x_minus_y(singlebolo_params,boloprops,savedir='./'):
    pixels = {}
    for b in boloprops:
        if b in singlebolo_params:
            pname = boloprops[b].physical_name
        else:
            continue
        pixel = pname.split('.')[0]
        if pixel not in pixels:
            pixels[pixel] = {}
            for band in BANDS:
                pixels[pixel][band] = []
        band = str(int(boloprops[b].band/core.G3Units.GHz))
        pixels[pixel][band].append(singlebolo_params[b]['rotPhi'])

    angdiff = {}
    for band in BANDS:
        angdiff[band] = []
    for pix in pixels:
        for band in BANDS:
            angs = pixels[pix][band]
            if len(angs) == 2:
                dang = CF.difference(angs[0],angs[1],deg=True,proj=2)
                angdiff[band].append(dang)
    pl.figure('Angle Diffs X-Y')
    for band in angdiff:
        pl.hist(angdiff[band],bins=np.linspace(0,180,61),color=BANDCOLORS[band],alpha=.5,label=band+' GHz: %.2f+-%.2f'%(np.mean(angdiff[band]),np.std(angdiff[band])))
    pl.grid()
    pl.title('X-Y bolo difference')
    pl.legend(fontsize=10)       
    pl.savefig(savedir+'X-Y_Best_Fit_AngleDiff_Hist.png')
    pl.close('all')


