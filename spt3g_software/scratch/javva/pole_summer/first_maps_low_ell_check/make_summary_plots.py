from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy, pylab
import matplotlib.pyplot as plt
import glob

plt.figure()
'''
Orginal files of PD and not, no map cuts
'''
a = pickle.load(open('ps_summary_LETLR_strict_cuts.pkl','rb')
a2 = pickle.load(open('ps_summary_notch_filter_strict_cuts.pkl','rb'))
a3 = pickle.load(open('ps_summary_p19_strict_cuts.pkl','rb'))


'''
a2 = pickle.load(open('ps_summary_LEPDLR_strict_cuts.pkl','rb'))
a3 = pickle.load(open('ps_summary_HETLR_strict_cuts.pkl','rb'))
a4 = pickle.load(open('ps_summary_HEPDLR_strict_cuts.pkl','rb'))
 
'''
k = 'T'
plt.figure()
for band in ['90','150','220']:
    plt.clf()
    for en, i in enumerate(a[k].keys()):
        
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko')
        plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko', label = 'Average of ell 46-298, low-ell T weights')
    plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*', label = 'average of ell 610-2456, low-ell T weights')

    for en, i in enumerate(a2[k].keys()):
        plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yo')
        plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*')

    plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yo', label = 'Average of ell 46-298, low-ell pair-diff weights')
    plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*', label = 'average of ell 610-2456, low-ell pair-diff weights')
    plt.legend()
    plt.semilogy()
    plt.xlabel('Days')
    plt.ylabel('Noise (uk-arcmin)')
    plt.title('Noise Averaging Down Over time - %s - %s'%(k,band))
    plt.ylim(55, 670500)
    plt.savefig('plots/%s_%s.png'%(band,k))
k = 'Q'
for band in ['90','150','220']:
    plt.clf()
    for en, i in enumerate(a[k].keys()):
        
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko')
        plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko', label = 'Average of ell 46-298, low-ell T weights')
    plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*', label = 'average of ell 610-2456, low-ell T weights')

    for en, i in enumerate(a2[k].keys()):
        plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yo')
        plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*')

    plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yo', label = 'Average of ell 46-298, low-ell pair-diff weights')
    plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*', label = 'average of ell 610-2456, low-ell pair-diff weights')
    plt.legend()
    plt.semilogy()
    plt.xlabel('Days')
    plt.ylabel('Noise (uk-arcmin)')
    plt.title('Noise Averaging Down Over time - %s - %s'%(k,band))
    plt.ylim(55, 10500)
    plt.savefig('plots/%s_%s.png'%(band,k))
k = 'U'
for band in ['90','150','220']:
    plt.clf()
    for en, i in enumerate(a[k].keys()):
        
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko')
        plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko', label = 'Average of ell 46-298, low-ell T weights')
    plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*', label = 'average of ell 610-2456, low-ell T weights')

    for en, i in enumerate(a2[k].keys()):
        plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yo')
        plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*')

    plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yo', label = 'Average of ell 46-298, low-ell pair-diff weights')
    plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*', label = 'average of ell 610-2456, low-ell pair-diff weights')
    plt.legend()
    plt.semilogy()
    plt.xlabel('Days')
    plt.ylabel('Noise (uk-arcmin)')
    plt.title('Noise Averaging Down Over time - %s - %s'%(k,band))
    plt.ylim(55, 10500)
    plt.savefig('plots/%s_%s.png'%(band,k))


'''
make summ plot
'''
plt.clf()
k = 'T'
for band in ['90','150','220']:
    plt.clf()
    k = 'T'
    for en, i in enumerate(a[k].keys()):
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ro')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ro', label = 'low-ell T weights, %s'%k)

    for en, i in enumerate(a2[k].keys()):
        plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'r*')
     #   plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*')

    plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'r*', label = 'low-ell pair-diff weights, %s'%k)
    #plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*', label = 'average of ell 610-2456, low-ell pair-diff weights')
    for en, i in enumerate(a3[k].keys()):
        plt.plot(en, np.mean(a3[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rD')
    plt.plot(en, np.mean(a3[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rD', label = 'high-ell T weights, %s'%k)
    for en, i in enumerate(a4[k].keys()):
        plt.plot(en, np.mean(a4[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rx')
    plt.plot(en, np.mean(a4[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rx', label = 'high-ell pair-diff weights, %s'%k)

    plt.legend()
  #  plt.semilogy()
    plt.xlabel('Days')
    plt.ylabel('Noise (uk-arcmin)')
    plt.title('Noise Averaging Down Over time - %s - %s'%(k,band))
    plt.ylim(55, 670500)
    #plt.savefig('plots/%s_%s.png'%(band,k))
    k = 'Q'
    #plt.clf()
    for en, i in enumerate(a[k].keys()):
        
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'go')
     #   plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'go', label = 'low-ell T weights, %s '%k)
    #plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*', label = 'average of ell 610-2456, low-ell T weights')

    for en, i in enumerate(a2[k].keys()):
        plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'g*')
     #   plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*')

    plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'g*', label = 'low-ell pair-diff weights, %s'%k)
    #plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*', label = 'average of ell 610-2456, low-ell pair-diff weights')
    for en, i in enumerate(a3[k].keys()):
        plt.plot(en, np.mean(a3[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gD')
    plt.plot(e, np.mean(a3[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gD', label = 'high-ell T weights, %s'%k)
    for en, i in enumerate(a4[k].keys()):
        plt.plot(en, np.mean(a4[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gx')
    plt.plot(en, np.mean(a4[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gx', label = 'high-ell pair-diff weights, %s'%k)


    plt.legend()
 #   plt.semilogy()
    plt.xlabel('Days')
    plt.ylabel('Noise (uk-arcmin)')
    plt.title('Noise Averaging Down Over time - %s - %s'%(k,band))
    plt.ylim(55, 10500)
    #plt.savefig('plots/%s_%s.png'%(band,k))
    k = 'U'
    #plt.clf()
    for en, i in enumerate(a[k].keys()):
        
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bo')
     #   plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bo', label = 'low-ell T weights, %s'%k)
    #plt.plot(en, np.mean(a[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k*', label = 'average of ell 610-2456, low-ell T weights')

    for en, i in enumerate(a2[k].keys()):
        plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'b*')
     #   plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*')

    plt.plot(en, np.mean(a2[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'b*', label = 'low-ell pair-diff weights, %s'%k)
    #plt.plot(en, np.mean(a2[k][i][band][50:])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y*', label = 'average of ell 610-2456, low-ell pair-diff weights')
    for en, i in enumerate(a3[k].keys()):
        plt.plot(en, np.mean(a3[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bD')
    plt.plot(en, np.mean(a3[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bD', label = 'high-ell T weights, %s'%k)
    for en, i in enumerate(a4[k].keys()):
        plt.plot(en, np.mean(a4[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bx')
    plt.plot(en, np.mean(a4[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bx', label = 'high-ell pair-diff weights, %s'%k)

    plt.legend( prop={'size': 6})
    plt.semilogy()
    plt.xlabel('Days')
    plt.ylabel('Noise (uk-arcmin)')
    plt.title('Avegage Noise between ell 46-298 vs Days Coadded,  %s GHz'%(band))
    plt.ylim(55, 670500)
    plt.savefig('plots/lowell_%s.png'%(band))



'''
next plot 
'''

plt.clf()
k = 'T'
for band in ['90','150','220']:
    plt.clf()
    k = 'T'
    for en, i in enumerate(a[k].keys()):
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ro')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ro', label = k)
    k = 'Q'
    for en, i in enumerate(a[k].keys()):
        
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'go')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'go', label = k)
    k = 'U'

    for en, i in enumerate(a[k].keys()):
        
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bo')
    plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bo', label = k)

    plt.legend()
    plt.semilogy()
    plt.xlabel('Days')
    plt.ylabel('Noise (uk-arcmin)')
    plt.title('Avegage Noise between ell 46-298 vs Days Coadded,  %s GHz'%(band))
    plt.ylim(55, 670500)
    plt.savefig('plots/let_lowell_%s.png'%(band))


'''
yet another plot
'''



band = '150'
for band in ['90','150','220']:
    #plt.clf()
    a = pickle.load(open('ADD_ps_summary_p19_strict_cuts.pkl','rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k', label ='TT')
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'r', label ='QQ')
    plt.plot(a['c'], (a['U'][list(a['U'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'g', label ='UU')
    a = pickle.load(open('ps_summary_all_maps_p19_strict_cuts_cf.pkl','rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k', label ='TT')
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'r', label ='QQ')
    plt.plot(a['c'], (a['U'][list(a['U'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'g', label ='UU')
    plt.semilogx()
    plt.semilogy()
    plt.legend()
    plt.xlabel('Ell')
    plt.ylim(55,500000)
    plt.title('Noise Power Spectrum T low-ell weights %s GHz'%band)
    plt.ylabel('Noise (uk-arcmin)')
#    plt.savefig('plots/%s_finalps.png'%band)

plt.plot(a['c'], (a2['T'][list(a2['T'].keys())[-1]]['150'])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k', label ='TT')
plt.plot(a['c'], (a2['Q'][list(a2['Q'].keys())[-1]]['150'])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'r', label ='QQ')
plt.plot(a['c'], (a2['U'][list(a2['U'].keys())[-1]]['150'])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'g', label ='UU')
plt.semilogx()
plt.semilogy()
plt.legend()
plt.xlabel('Ell')
plt.ylabel('Noise (uk-arcmin)')


t,q,u = mapmaker.mapmakerutils.remove_weight(i['T'], i['Q'], i['U'],i['
    ...: Wpol'\                                
    ...: ]) 
apod = mapspectra.apodmask.make_border_apodization((i['Wpol']),apod_type='cos', radius_arcmin=90.)
plt.figure(3)
plt.clf()
plt.imshow(t*apod)
plt.title('T - 30 Day Coadd')
plt.colorbar()
plt.clim(-1,1)
plt.figure(4)
plt.clf()
plt.imshow(q*apod)
plt.title('Q - 30 Day Coadd')
plt.colorbar()
plt.clim(-.1,.1)
plt.figure(5)
plt.clf()
plt.imshow(u*apod)
plt.title('U - 30 Day Coadd')
plt.colorbar()
plt.clim(-.1,.1)


plt.figure(6)
plt.plot(np.asarray(i['Wpol'].TT)[:,1147], label = 'T')
plt.plot(np.asarray(i['Wpol'].QQ)[:,1147], label = 'Q')
plt.plot(np.asarray(i['Wpol'].UU)[:,1147], label = 'U')


'''
Make Overplotted Power Spectrum Plot
'''

band = '150'
for band in ['90','150','220']:
    plt.clf()
    nm = 'p10'
    a = pickle.load(open('ps_summary_all_maps_%s_strict_cuts_cf.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k--', label ='noise TT - %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k--', label ='noise QQ - %s '%nm)
    a = pickle.load(open('ADD_ps_summary_all_maps_%s_strict_cuts.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k', label ='TT %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'k', label ='QQ %s'%nm)
    nm = 'p19'
    a = pickle.load(open('ps_summary_all_maps_%s_strict_cuts_cf.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m--', label ='noise TT - %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m--', label ='noise QQ - %s '%nm)
    a = pickle.load(open('ADD_ps_summary_all_maps_%s_strict_cuts.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m', label ='TT %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m', label ='QQ %s'%nm)
    nm = 'p4'
    a = pickle.load(open('ps_summary_LETLR_strict_cuts.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band]/2)**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y--', label ='noise TT - %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band]/2)**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y--', label ='noise QQ - %s '%nm)
    a = pickle.load(open('ADD_ps_summary_LETLR_strict_cuts_strict_cuts.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y', label ='TT %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'y', label ='QQ %s'%nm)


plt.clf()
k = 'Q'
for en, i in enumerate(a[k].keys()):
    en = i.split('_')[-1].split('.g3')[0]
    if '44' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yx')
    if '52' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rx')
    if '67' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bx')
    if '59' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gx')
plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'kx', label= 'Q')

k = 'T'
for en, i in enumerate(a[k].keys()):
    en = i.split('_')[-1].split('.g3')[0]
    if '44' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'yo')
    if '52' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ro')
    if '67' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bo')
    if '59' in i.split('/')[5]:
        plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'go')
plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ko', label = 'T')
plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'ys', label = '44.75')
plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rs', label = '52.25')
plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bs', label = '67.25')
plt.plot(en, np.mean(a[k][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gs', label = '59.75')

plt.legend()
                plt.title('Noise for Individual Observations in T and Q - %s GHz'%band)
plt.xlabel('Observation Number')
plt.semilogy()
plt.ylim(10e1,10e5)
plt.ylabel('Mean Noise of Ell = 46-430 (uK-arcmin)')
plt.savefig('plots/wiki/%s_indyfields.png'%band)


field = '44.75'
plt.clf()
for field in ['44.75','52.25','67.25','59.75']:
    plt.clf()
    for en, i in enumerate(a[k].keys()): 
        if field in i.split('/')[5]:
            band = '90'
            plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rx')
            band = '150'
            plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bx')
            band = '220'
            plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gx')
    plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rx', label = '90')
    plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bx', label = '150')
    plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gx', label = '220')
    plt.legend()
    plt.semilogx()
    plt.semilogy()
    plt.xlim(10e2,10e5)
    plt.ylim(10e1,10e5)
    plt.xlabel('T Noise (uK-arcmin)')
    plt.ylabel('Q Noise(uK-arcmin)')
    plt.title('T to P correlation for %s Subfield'%field)
    plt.savefig('plots/wiki/sfcorr%s.png'%field)

for field in ['44.75','52.25','67.25','59.75']:
    plt.clf()
    for en, i in enumerate(a[k].keys()): 
        if field in i.split('/')[5]:
            band = '90'
            plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rx')
            band = '150'
            plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bx')
            band = '220'
            plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gx')
    plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'rx', label = '90')
    plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'bx', label = '150')
    plt.plot(np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK), np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'gx', label = '220')
    plt.legend()
    plt.semilogx()
    plt.semilogy()
    plt.xlim(10e2,10e5)
    plt.ylim(10e1,10e5)
    plt.xlabel('T Noise (uK-arcmin)')
    plt.ylabel('Q Noise(uK-arcmin)')
    plt.title('T to P correlation for %s Subfield'%field)
    plt.savefig('plots/wiki/sfcorr%s.png'%field)


for field in ['44.75','52.25','67.25','59.75']:
    plt.clf()
    for l in [1]:
                band = '90'
                hsts = []
                for en, i in enumerate(a[k].keys()):
                    if field in i.split('/')[5]:
                        hsts.append((np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/( np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK)))
                plt.hist(hsts, histtype = 'step', color = 'r',bins = np.linspace(0,.1,50), label = band)
                band = '150'
                hsts = []
                for en, i in enumerate(a[k].keys()):
                    if field in i.split('/')[5]:
                        hsts.append((np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/( np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK)))
                plt.hist(hsts, histtype = 'step', color = 'b',bins = np.linspace(0,.1,50), label = band)
                band = '220'
                hsts = []
                for en, i in enumerate(a[k].keys()):
                    if field in i.split('/')[5]:
                        hsts.append((np.mean(a['Q'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK))/( np.mean(a['T'][i][band][3:35])**0.5 / (core.G3Units.arcmin * core.G3Units.uK)))
                plt.hist(hsts, histtype = 'step', color = 'g',bins = np.linspace(0,.1,50), label = band)
 
    plt.legend()
    plt.xlabel('Q Noise/ T Noise (proxy for fractional T leakage)')
    plt.ylabel('Number of Observations')
    plt.title('T -> P leakage for %s Field'%field)
    plt.savefig('plots/wiki/leakage_%s.png'%field)


    nm = 'p19'
    a = pickle.load(open('ps_summary_all_maps_%s_strict_cuts_cf.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m--', label ='noise TT - %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m--', label ='noise QQ - %s '%nm)
    a = pickle.load(open('ADD_ps_summary_%s_strict_cuts.pkl'%nm,'rb'))
    plt.plot(a['c'], (a['T'][list(a['T'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m', label ='TT %s'%nm)
    plt.plot(a['c'], (a['Q'][list(a['Q'].keys())[-1]][band])**0.5 / (core.G3Units.arcmin * core.G3Units.uK),'m', label ='QQ %s'%nm)
