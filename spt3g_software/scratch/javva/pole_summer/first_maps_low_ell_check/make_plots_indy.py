from spt3g import core, coordinateutils, mapmaker, mapspectra
import numpy, pickle, sys, copy, pylab
from spt3g.mapspectra import map_analysis
import matplotlib.pylab as plt
import glob

for band in ['90','150','220']:
        band = s.split('_')[2]
        a = pickle.load(open(s,'rb'))
        plt.clf()
        for field in ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']:
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field:
                                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'ro')
                for i in a['T'].keys():
                        if i.split('/')[-4] ==field:
                                plt.plot(i.split('/')[-2], (a['T'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'go')
                for i in a['U'].keys():
                        if i.split('/')[-4] ==field:
                                plt.plot(i.split('/')[-2], (a['U'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'bo')
                plt.semilogy()
                plt.xlabel('Observation Number')
                plt.ylabel('Noise at ell = 100 (uK-arcmin)')
                plt.title(field + ' '+ band)
                plt.plot(i.split('/')[-2], (a['U'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'go', label = 'T')
                plt.plot(i.split('/')[-2], (a['U'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'ro', label = 'Q')
                plt.plot(i.split('/')[-2], (a['U'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'bo',label = 'U')
                plt.savefig('plots/summaryplots/%s_%s.png'%(field,band))
                plt.xlim(7e7, 7.4e7)
                plt.ylim(1e2,1e7)
                plt.savefig('plots/summaryplots/zoom_%s_%s.png'%(field,band))


for s in ['ps_summary_90_landr_newana.pkl','ps_summary_150_landr_newana.pkl', 'ps_summary_220_landr_newana.pkl']:
        band = s.split('_')[2]
        a = pickle.load(open(s,'rb'))
        plt.clf()
        for i in [1]:
                field = 'ra0hdec-44.75'
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'ro')
                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),'ro', label = '%s'%field)
                field = 'ra0hdec-52.25'
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'y*')
                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),'y*', label = '%s'%field)
                field = 'ra0hdec-59.75'
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'kD')
                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),'kD', label = '%s'%field)
                field = 'ra0hdec-67.25'
                
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),'gx')
                plt.plot(i.split('/')[-2], (a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),'gx', label = '%s'%field)
                plt.semilogy()
                plt.xlabel('Observation Number')
                plt.ylabel('Noise at ell = 100 (uK-arcmin)')
                plt.title('Field Comparison Q '+ band)
                plt.savefig('plots/summaryplots/fc_Q_%s_just7.png'%(band))
                plt.xlim(7e7, 7.4e7)
                plt.ylim(1e2,1e7)
                plt.savefig('plots/summaryplots/zoom_fc_Q_%s_just7.png'%(band))


for s in ['ps_summary_90_landr_newana.pkl','ps_summary_150_landr_newana.pkl', 'ps_summary_220_landr_newana.pkl']:
        band = s.split('_')[2]
        a = pickle.load(open(s,'rb'))
        plt.clf()
        for i in [1]:
                field = 'ra0hdec-44.75'
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'ro')
                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'ro', label = '%s'%field)
                field = 'ra0hdec-52.25'
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'y*')
                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'y*', label = '%s'%field)

                field = 'ra0hdec-59.75'
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'kD')
                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'kD', label = '%s'%field)

                field = 'ra0hdec-67.25'
                for i in a['Q'].keys():
                        if i.split('/')[-4] ==field and i.split('/')[-2][0] == '7':
                                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin * core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'gx')
                plt.plot((a['Q'][i][8])**.5/(core.G3Units.arcmin *core.G3Units.uK),(np.mean(a['Q'][i][50:]))**.5/(core.G3Units.arcmin * core.G3Units.uK),'gx', label = '%s'%field)
                

                plt.semilogy()
                plt.semilogy()
                plt.ylabel('Average value between ell = 610 - 2500 (uK-arcmin)')
                plt.xlabel('Noise at ell = 100 (uK-arcmin)')
                plt.title('Does Low Ell correllate to High ell -- '+ band)
                plt.savefig('plots/summaryplots/fc_Q_%s_just7.png'%(band))
                plt.xlim(0, 10000)
                plt.ylim(0,10000)
                plt.savefig('plots/summaryplots/zoom_fc_Q_%s_just7.png'%(band))
