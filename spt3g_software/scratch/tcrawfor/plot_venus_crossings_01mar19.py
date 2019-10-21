import pickle
import matplotlib.backends.backend_pdf

#outdict = pickle.load(open('/spt/user/tcrawfor/public/venus_fit_out_28feb19.pkl'))

#gnames=[name for name in names if outdict[name]['left']['pars_full'][0] > 200. and outdict[name]['left']['pars_full'][2] < 5.*0.42 and outdict[name]['left']['pars_lead'][2] > 0.8*0.42 and outdict[name]['left']['pars_trail'][2] > 0.8*0.42 and outdict[name]['left']['pars_full'][2] > 0.8*0.42]

pnames = [bp[name].physical_name for name in gnames]
spnames = np.argsort(np.asarray(pnames))

pdf = matplotlib.backends.backend_pdf.PdfPages('/spt/user/tcrawfor/public/venus_crossings_01mar19.pdf')

#figure(figsize=(15,12))
npx = 3
npy = 3
nptot = np.int(npx*npy)

el_venus = 20.*np.pi/180.
scanspeed_onsky = 0.3*np.cos(el_venus)
arcmin_per_sample = 60.*scanspeed_onsky/152.6
nx = 300
xtemp = (np.arange(nx)-nx/2)*arcmin_per_sample

for q in np.arange(len(gnames)):
    qq = spnames[q]
    gnq = gnames[qq]
    pnq = pnames[qq]
    odnl = outdict[gnq]['left']
    pnum = q - (q/nptot)*nptot
    subplot(npx, npy, pnum+1)
    plot(xtemp,odnl['data'],'o')
    xlim(-5,5)
    if pnum >= nptot - npx:
        xlabel('offset from source center [arcmin]')
    title(pnq)
    plot(xtemp,odnl['fit_full'],label='FWHM='+"{:4.2f}".format(odnl['pars_full'][2]/0.42))
    plot(xtemp,odnl['fit_lead'],label='FWHM='+"{:4.2f}".format(odnl['pars_lead'][2]/0.42))
    plot(xtemp,odnl['fit_trail'],label='FWHM='+"{:4.2f}".format(odnl['pars_trail'][2]/0.42))
    legend(loc=1,fontsize='x-small')
    if pnum == nptot-1:
        pdf.savefig()
        clf()

pdf.close()
