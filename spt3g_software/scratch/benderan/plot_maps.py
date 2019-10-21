from spt3g import core, mapmaker, coordinateutils
import matplotlib.pyplot as plt
import numpy as np
import os,sys
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import pdb
from scipy import ndimage
import time

def plot_map_bands(mapfile,xlm=[],ylm=[],tickoff=True,colorbar_on=False,plotbeam=False,title_inpanel=True, anchor_clim=0.,unpol=True):
    data=list(core.G3File(mapfile))
    if unpol:
        unweighted90=data[0]['T']/data[0]['Wunpol'].TT
        unweighted150=data[1]['T']/data[1]['Wunpol'].TT
        unweighted220=data[2]['T']/data[2]['Wunpol'].TT
    else:
        unweighted90=data[0]['T']/data[0]['Wpol'].TT
        unweighted150=data[2]['T']/data[1]['Wpol'].TT
        unweighted220=data[1]['T']/data[2]['Wpol'].TT
    #clims=[np.max(unweighted90),np.max(unweighted150),np.max(unweighted220)]

    mapres=data[0]['T'].res/core.G3Units.arcmin    

    if len(xlm)==0:
        xlm=[0,np.shape(unweighted90)[0]]
    if len(ylm)==0:
        ylm=[0,np.shape(unweighted90)[0]]


    circ_cent=((xlm[1]-xlm[0])*0.125+xlm[0],(ylm[1]-ylm[0])*0.125+ylm[0])
    circ_rad=0.5/mapres
    print circ_cent, circ_rad
    circ1=plt.Circle(circ_cent,circ_rad,color='w',fill=False,linewidth=1.5)
    

    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)
    pp=plt.imshow(unweighted90,interpolation='none')
    plt.xlim(xlm)
    plt.ylim(ylm)
    if tickoff:
        xtk=[]
        ytk=[]
    else:
        xtk=plt.xticks()[0]
        ytk=plt.yticks()[0]
    plt.xticks(xtk)
    plt.yticks(ytk)
    if title_inpanel:
        plt.text((xlm[1]-xlm[0])*0.7+xlm[0], (ylm[1]-ylm[0])*0.875+ylm[0],'90 GHz', color='w', fontsize=15,fontweight='semibold')
    clm=pp.get_clim()
    print clm
    if anchor_clim !=0.:
        plt.clim([anchor_clim,clm[1]])
    plt.title('90 GHz')
    if colorbar_on:
        plt.colorbar()
    plt.subplot(1,3,2)
    pp=plt.imshow(unweighted150,interpolation='none')
    plt.xlim(xlm)
    plt.ylim(ylm)
    plt.xticks(xtk)
    plt.yticks(ytk)
    clm=pp.get_clim()
    print clm
    if anchor_clim !=0.:
        plt.clim([anchor_clim,clm[1]])
    ax=plt.gca()
    ax.add_artist(circ1)
    plt.text(circ_cent[0]+circ_rad*2, circ_cent[1]-circ_rad*0.8,'1\'',color='w',fontsize=15,fontweight='medium') 
    if title_inpanel:
        plt.text((xlm[1]-xlm[0])*0.65+xlm[0], (ylm[1]-ylm[0])*0.875+ylm[0],'150 GHz', color='w', fontsize=15,fontweight='semibold')
    #plt.clim([0,np.max(clims)])
    plt.title('150 GHz')
    if colorbar_on:
        plt.colorbar()
    plt.subplot(1,3,3)
    pp=plt.imshow(unweighted220,interpolation='none')
    plt.xlim(xlm)
    plt.ylim(ylm)
    plt.xticks(xtk)
    plt.yticks(ytk)
    clm=pp.get_clim()
    print clm
    if anchor_clim !=0.:
        plt.clim([anchor_clim,clm[1]])
    #plt.clim([0,np.max(clims)])
    if title_inpanel:
        plt.text((xlm[1]-xlm[0])*0.65+xlm[0], (ylm[1]-ylm[0])*0.875+ylm[0],'220 GHz', color='w', fontsize=15,fontweight='semibold')
    plt.title('220 GHz')
    if colorbar_on:
        plt.colorbar()
    return


def plot_map_bands_logscale(mapfile):
    data=list(core.G3File(mapfile))
    unweighted90=data[0]['T']/data[0]['Wunpol'].TT
    unweighted150=data[1]['T']/data[1]['Wunpol'].TT
    unweighted220=data[2]['T']/data[2]['Wunpol'].TT
    plt.figure()
    plt.subplot(1,3,1)
    plt.imshow(np.log10(unweighted90),interpolation='none')
    plt.title('90 GHz')
    plt.colorbar()
    plt.subplot(1,3,2)
    plt.imshow(np.log10(unweighted150),interpolation='none')
    plt.title('150 GHz')
    plt.colorbar()
    plt.subplot(1,3,3)
    plt.imshow(np.log10(unweighted220),interpolation='none')
    plt.title('220 GHz')
    plt.colorbar()
    return


#mapfile='/spt/user/ddutcher/coadds/2018_500d_coadd.g3'
mapfile='/spt/user/ddutcher/coadds/w1-4Hz_20190310.g3'
def plot_500d_field_map(mapfile,xlm=[],ylm=[]):
    data=list(core.G3File(mapfile))
    unweighted90=data[1]['T']/data[1]['Wpol'].TT
    unweighted150=data[2]['T']/data[2]['Wpol'].TT
    unweighted220=data[3]['T']/data[3]['Wpol'].TT
    #clims=[np.max(unweighted90),np.max(unweighted150),np.max(unweighted220)]                                                       

    mapres=data[1]['T'].res/core.G3Units.arcmin

    #if len(xlm)==0:
    #    xlm=[0,np.shape(unweighted90)[0]]
    #if len(ylm)==0:
    #    ylm=[0,np.shape(unweighted90)[0]]

    tmp150 = np.asarray(unweighted150)
    tmp150weight=np.asarray(data[2]['Wpol'].TT)
    apod_mask=np.zeros(tmp150.shape)
    apod_mask[np.where(tmp150weight > 80)]=1.
    tmp150=tmp150*apod_mask
    tmp150[np.isnan(tmp150)] = 0
    tmp150[np.where(tmp150 ==0)]=-100
    
    #plt.figure(figsize=(14,7))
    plt.figure()
    plt.imshow(tmp150,interpolation='none',cmap='bone')
    plt.clim(-8,8)
    plt.xticks([])
    plt.yticks([])
#    plt.ylim(0,1280)
    plt.text(250,850,'150 GHz', fontsize=24,color='w',weight='bold')
    plt.savefig('field_500d_150ghz_black.png')

    tmp90 = np.asarray(unweighted90)
    tmp90weight=np.asarray(data[1]['Wpol'].TT)
    #apod_mask=np.zeros(tmp90.shape)
    #apod_mask[np.where(tmp90weight > 80)]=1.
    tmp90=tmp90*apod_mask
    tmp90[np.isnan(tmp90)] = 0
    tmp90[np.where(tmp90 ==0)]=-100

    plt.figure(figsize=(14,7))
    plt.imshow(tmp90,interpolation='none',cmap='bone')
    plt.clim(-0.5,0.5)
    plt.xticks([])
    plt.yticks([])
    plt.text(250,850,'95 GHz', fontsize=24,color='w',weight='bold')
 #   plt.ylim(0,1280)
    plt.savefig('field_500d_90ghz_black.png')

    '''
    plt.figure(figsize=(14,7))
    plt.imshow(tmp90,interpolation='none',cmap='bone')
    plt.clim(-0.5,0.5)
    plt.xticks([])
    plt.yticks([])
    plt.text(250,850,'95 GHz', fontsize=24,color='w',weight='bold')
    plt.ylim(0,1280)
    rectangle = plt.Rectangle((850, 570), 500, 270, fc='none',ec='w',linewidth=3)
    plt.gca().add_patch(rectangle)
    plt.savefig('field_500d_90ghz_black_rectangle.png')

    plt.figure(figsize=(14,7))
    plt.imshow(tmp90,interpolation='none',cmap='bone')
    plt.clim(-0.5,0.5)
    plt.xticks([])
    plt.yticks([])
    plt.ylim(570,840)
    plt.xlim(850,1350)
    plt.savefig('field_500d_90ghz_black_subfield.png')
    '''

def plot_field_map(mapfile,saveprefix='field_1500d_'):
    data=list(core.G3File(mapfile))
    unweighted90=data[1]['T']/data[1]['Wpol'].TT
    unweighted150=data[2]['T']/data[2]['Wpol'].TT
    unweighted220=data[3]['T']/data[3]['Wpol'].TT
    mapres=data[1]['T'].res/core.G32zUnits.arcmin

    ra,dec=get_map_ra_dec(data[1]['T'])
    ra_contour_levels=np.arange(-5,6)
    dec_contour_levels=np.arange(-70,-20,10)

    tmp150 = np.asarray(unweighted150)
    tmp150weight=np.asarray(data[2]['Wpol'].TT)
    apod_mask=np.zeros(tmp150.shape)
    apod_mask[np.where(tmp150weight > 80)]=1.
    tmp150=tmp150*apod_mask
    tmp150[np.isnan(tmp150)] = 0
    tmp150[np.where(tmp150 ==0)]=-100

    plt.figure(figsize=(14,7))
    plt.imshow(tmp150,interpolation='none',cmap='bone')
    plt.clim(-.2,.2)

    ra_vec=ra[2700,:]
    ra_tck_pos=[]
    ra_tck_levels=[]
    for kk in ra_contour_levels:
        tck_tmp=ra_vec[np.where( np.abs(ra_vec-kk)==np.min(np.abs(ra_vec-kk)))[0][0]]
        if np.abs(tck_tmp-kk) < 0.1:
            ra_tck_pos.append(np.where( np.abs(ra_vec-kk)==np.min(np.abs(ra_vec-kk)))[0][0])
            if tck_tmp < 0:
                tck_tmp=tck_tmp+24
            ra_tck_levels.append(str(np.int32(np.round(tck_tmp)))+'$\,$h')
    
    dec_vec=dec[:,50]
    dec_vec2=dec[:,60]
    dec_tck_pos=[]
    dec_tck_levels=[]
    dec_tck_angles=[]
    for kk in dec_contour_levels:
        tck_tmp=dec_vec[np.where( np.abs(dec_vec-kk)==np.min(np.abs(dec_vec-kk)))[0][0]]
        if np.abs(tck_tmp-kk) < 0.1:
            dec_tck_pos.append(np.where( np.abs(dec_vec-kk)==np.min(np.abs(dec_vec-kk)))[0][0])
            dec_tck_levels.append(str(np.int32(np.round(tck_tmp)))+'$^{\circ}$')
    for ii,kk in enumerate(dec_contour_levels):
        tck_tmp=np.where( np.abs(dec_vec2-kk)==np.min(np.abs(dec_vec2-kk)))[0][0]
        if np.abs(tck_tmp-kk) < 0.1:
            dec_tck_angles.append(np.asin(10./(tck_tmp-dec_tck_pos[ii]))*180./np.pi)
    print dec_tck_angles

    plt.yticks(dec_tck_pos,dec_tck_levels,color='w',fontsize=12)
    plt.ylim(100,2800)
    plt.xticks([])
    plt.yticks([])
    for ii,kk in enumerate(ra_tck_pos):
        plt.text(kk+50,2700,ra_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'))
    for ii,kk in enumerate(dec_tck_pos):
        plt.text(40,kk+200-ii*20,dec_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'),rotation=45-ii*6)
    plt.contour(ra,ra_contour_levels,colors='w', linestyles='dashed',linewidths=.75)
    plt.contour(dec,dec_contour_levels,colors='w',linestyles='dashed',linewidths=.75)
    #plt.text(250,2500,'150 GHz', fontsize=24,color='w',weight='bold',bbox=dict(facecolor='black'))
    plt.text(3600,250,'150 GHz', fontsize=24,color='w',weight='bold',bbox=dict(facecolor='black'))
    plt.savefig(saveprefix+'150ghz_black.png')

    tmp90 = np.asarray(unweighted90)
    tmp90weight=np.asarray(data[1]['Wpol'].TT)
    #apod_mask=np.zeros(tmp90.shape)                                                                       
    #apod_mask[np.where(tmp90weight > 80)]=1.                                                              
    tmp90=tmp90*apod_mask
    tmp90[np.isnan(tmp90)] = 0
    tmp90[np.where(tmp90 ==0)]=-100


    #### 95 Ghz TT
    plt.figure(figsize=(14,7))
    plt.imshow(tmp90,interpolation='none',cmap='bone')
    plt.clim(-0.2,0.2)
    plt.ylim(100,2800)

    plt.ylim(100,2800)
    plt.xticks([])
    plt.yticks([])
    for ii,kk in enumerate(ra_tck_pos):
        plt.text(kk+50,2700,ra_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'))
    for ii,kk in enumerate(dec_tck_pos):
        plt.text(40,kk+200-ii*20,dec_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'),rotation=45-ii*6)
    plt.contour(ra,ra_contour_levels,colors='w', linestyles='dashed',linewidths=1)
    plt.contour(dec,dec_contour_levels,colors='w',linestyles='dashed',linewidths=1)
    plt.text(3600,250,'95 GHz', fontsize=24,color='w',weight='bold',bbox=dict(facecolor='black'))
        
    plt.savefig(saveprefix+'95ghz_black.png')

    ## 220 Ghz TT
    
    tmp220 = np.asarray(unweighted220)
    tmp220weight=np.asarray(data[1]['Wpol'].TT)
    
    tmp220=tmp220*apod_mask
    tmp220[np.isnan(tmp220)] = 0
    tmp220[np.where(tmp220 ==0)]=-100


    plt.figure(figsize=(14,7))
    plt.imshow(tmp220,interpolation='none',cmap='bone')
    plt.clim(-0.3,0.3)
    plt.ylim(100,2800)

    plt.ylim(100,2800)
    plt.xticks([])
    plt.yticks([])
    for ii,kk in enumerate(ra_tck_pos):
        plt.text(kk+50,2700,ra_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'))
    for ii,kk in enumerate(dec_tck_pos):
        plt.text(40,kk+200-ii*20,dec_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'),rotation=45-ii*6)
    plt.contour(ra,ra_contour_levels,colors='w', linestyles='dashed',linewidths=1)
    plt.contour(dec,dec_contour_levels,colors='w',linestyles='dashed',linewidths=1)
    plt.text(3600,250,'220 GHz', fontsize=24,color='w',weight='bold',bbox=dict(facecolor='black'))

    plt.savefig(saveprefix+'220ghz_black.png')
    


    return

    

def get_map_ra_dec(map):
    #map=data[1]['T']
    ra=np.zeros(np.transpose(map.shape))
    dec=np.zeros(np.transpose(map.shape))
    
    for kk in range(map.shape[0]):
        for jj in range(map.shape[1]):
            ra[kk,jj]=map.pixel_to_angle(jj,kk)[0]/core.G3Units.rahr
            dec[kk,jj]=map.pixel_to_angle(jj,kk)[1]/core.G3Units.deg
    return ra, dec



def plot_field_map_general(mapfile,map_index=1,map_type='', label='95 GHz',clim=[-8,8],savedir='/home/benderan',saveprefix='field_1500d_95GHz',smooth=False,plot_coords=True,bckgrnd='black'):
    data=list(core.G3File(mapfile))
    if map_type=='T':
        weight_map=data[map_index]['Wpol'].TT
    elif map_type=='Q':
        weight_map=data[map_index]['Wpol'].QQ
    elif map_type=='U':
        weight_map=data[map_index]['Wpol'].UU

    if map_type=='T' and smooth==True:
        raw_input('are you sure you want to smooth the T map?')

    #    if map_type=='T':
    unweighted=data[map_index][map_type]/weight_map
    #    else:
    #        unweighted=data[map_index][map_type]
    print('done reading in the map file')
    mapres=data[map_index][map_type].res/core.G3Units.arcmin

    ra,dec=get_map_ra_dec(data[map_index][map_type])
    ra_contour_levels=np.arange(-5,6)
    dec_contour_levels=np.arange(-70,-20,10)

    if smooth:
        unweighted=ndimage.gaussian_filter(unweighted,4)

    tmp150 = np.asarray(unweighted)
    tmp150weight=np.asarray(weight_map)
    apod_mask=np.zeros(tmp150.shape)
    apod_mask[np.where(tmp150weight > 80)]=1.
    tmp150=tmp150*apod_mask
    tmp150[np.isnan(tmp150)] = 0
    if bckgrnd=='black':
        tmp150[np.where(tmp150 ==0)]=-100
        contcol='w'
    elif bckgrnd=='white':
        tmp150[np.where(tmp150 ==0)]=100
        contcol='k'

    #ra=ra*(1.-apod_mask)
    #dec=dec*(1.-apod_mask)
    plt.figure(figsize=(14,7))
    ax=plt.imshow(tmp150,interpolation='none',cmap='bone')
    time.sleep(1)
    ax.set_clim(clim)
    

    ra_vec=ra[2700,:]
    ra_tck_pos=[]
    ra_tck_levels=[]
    for kk in ra_contour_levels:
        tck_tmp=ra_vec[np.where( np.abs(ra_vec-kk)==np.min(np.abs(ra_vec-kk)))[0][0]]
        if np.abs(tck_tmp-kk) < 0.1:
            ra_tck_pos.append(np.where( np.abs(ra_vec-kk)==np.min(np.abs(ra_vec-kk)))[0][0])
            if tck_tmp < 0:
                tck_tmp=tck_tmp+24
            ra_tck_levels.append(str(np.int32(np.round(tck_tmp)))+'$\,$h')

    dec_vec=dec[:,50]
    dec_vec2=dec[:,60]
    dec_tck_pos=[]
    dec_tck_levels=[]
    dec_tck_angles=[]
    for kk in dec_contour_levels:
        tck_tmp=dec_vec[np.where( np.abs(dec_vec-kk)==np.min(np.abs(dec_vec-kk)))[0][0]]
        if np.abs(tck_tmp-kk) < 0.1:
            dec_tck_pos.append(np.where( np.abs(dec_vec-kk)==np.min(np.abs(dec_vec-kk)))[0][0])
            dec_tck_levels.append(str(np.int32(np.round(tck_tmp)))+'$^{\circ}$')
    for ii,kk in enumerate(dec_contour_levels):
        tck_tmp=np.where( np.abs(dec_vec2-kk)==np.min(np.abs(dec_vec2-kk)))[0][0]
        if np.abs(tck_tmp-kk) < 0.1:
            dec_tck_angles.append(np.asin(10./(tck_tmp-dec_tck_pos[ii]))*180./np.pi)
    print(dec_tck_angles)

    plt.yticks(dec_tck_pos,dec_tck_levels,color=contcol,fontsize=12)
    plt.ylim(100,2800)
    plt.xticks([])
    plt.yticks([])
    if plot_coords:
        for ii,kk in enumerate(ra_tck_pos):
                plt.text(kk+60,2670,ra_tck_levels[ii],color=contcol,fontsize=12,bbox=dict(facecolor=bckgrnd,edgecolor=bckgrnd))
        for ii,kk in enumerate(dec_tck_pos):
            plt.text(40,kk+200-ii*20,dec_tck_levels[ii],color=contcol,fontsize=12,bbox=dict(facecolor=bckgrnd,edgecolor=bckgrnd),rotation=45-ii*6)
        plt.contour(ra,ra_contour_levels,colors=contcol, linestyles='dashed',linewidths=.75)
        plt.contour(dec,dec_contour_levels,colors=contcol,linestyles='dashed',linewidths=.75)
        
    plt.text(3600,250,label, fontsize=24,color=contcol,weight='bold',bbox=dict(facecolor=bckgrnd,edgecolor=bckgrnd))
    plt.tight_layout()
    plt.savefig(os.path.join(savedir,saveprefix+'_'+bckgrnd+'.png'))
    return 

def plot_map_and_zoom(mapfile,map_index=1,map_type='', label='95 GHz',clim=[-8,8],savedir='/home/benderan',saveprefix='field_1500d_95GHz',smooth=False,plot_coords=True,zoom_bounds=[]):
    '''
    zoom_bounds: [[x1,x2],[y1,y2]]
    
    '''

    data=list(core.G3File(mapfile))
    if map_type=='T':
        weight_map=data[map_index]['Wpol'].TT
    elif map_type=='Q':
        weight_map=data[map_index]['Wpol'].QQ
    elif map_type=='U':
        weight_map=data[map_index]['Wpol'].UU

    if map_type=='T' and smooth==True:
        raw_input('are you sure you want to smooth the T map?')

    unweighted=data[map_index][map_type]/weight_map
    mapres=data[map_index][map_type].res/core.G3Units.arcmin

    ra,dec=get_map_ra_dec(data[map_index][map_type])
    ra_contour_levels=np.arange(-5,6)
    dec_contour_levels=np.arange(-70,-20,10)

    if smooth:
        unweighted=ndimage.gaussian_filter(unweighted,4)

    tmp150 = np.asarray(unweighted)
    tmp150weight=np.asarray(weight_map)
    apod_mask=np.zeros(tmp150.shape)
    apod_mask[np.where(tmp150weight > 80)]=1.
    tmp150=tmp150*apod_mask
    tmp150[np.isnan(tmp150)] = 0
    tmp150[np.where(tmp150 ==0)]=-100
    plt.figure(figsize=(14,7))
    ax=plt.imshow(tmp150,interpolation='none',cmap='bone')
    ax.set_clim(clim)


    ra_vec=ra[2700,:]
    ra_tck_pos=[]
    ra_tck_levels=[]
    for kk in ra_contour_levels:
        tck_tmp=ra_vec[np.where( np.abs(ra_vec-kk)==np.min(np.abs(ra_vec-kk)))[0][0]]
        if np.abs(tck_tmp-kk) < 0.1:
            ra_tck_pos.append(np.where( np.abs(ra_vec-kk)==np.min(np.abs(ra_vec-kk)))[0][0])
            if tck_tmp < 0:
                tck_tmp=tck_tmp+24
            ra_tck_levels.append(str(np.int32(np.round(tck_tmp)))+'$\,$h')

    dec_vec=dec[:,50]
    dec_vec2=dec[:,60]
    dec_tck_pos=[]
    dec_tck_levels=[]
    dec_tck_angles=[]
    for kk in dec_contour_levels:
        tck_tmp=dec_vec[np.where( np.abs(dec_vec-kk)==np.min(np.abs(dec_vec-kk)))[0][0]]
        if np.abs(tck_tmp-kk) < 0.1:
            dec_tck_pos.append(np.where( np.abs(dec_vec-kk)==np.min(np.abs(dec_vec-kk)))[0][0])
            dec_tck_levels.append(str(np.int32(np.round(tck_tmp)))+'$^{\circ}$')
    for ii,kk in enumerate(dec_contour_levels):
        tck_tmp=np.where( np.abs(dec_vec2-kk)==np.min(np.abs(dec_vec2-kk)))[0][0]
        if np.abs(tck_tmp-kk) < 0.1:
            dec_tck_angles.append(np.asin(10./(tck_tmp-dec_tck_pos[ii]))*180./np.pi)
    print dec_tck_angles

    plt.yticks(dec_tck_pos,dec_tck_levels,color='w',fontsize=12)
    plt.ylim(100,2800)
    plt.xticks([])
    plt.yticks([])
    if plot_coords:
        for ii,kk in enumerate(ra_tck_pos):
            plt.text(kk+50,2700,ra_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'))
        for ii,kk in enumerate(dec_tck_pos):
            plt.text(40,kk+200-ii*20,dec_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'),rotation=45-ii*6)
        plt.contour(ra,ra_contour_levels,colors=bckgnd, linestyles='dashed',linewidths=.75)
        plt.contour(dec,dec_contour_levels,colors=bckgnd,linestyles='dashed',linewidths=.75)

    plt.text(3600,250,label, fontsize=24,color=bckgnd,weight='bold',bbox=dict(facecolor='black'))
    rect=patches.Rectangle([zoom_bounds[0][0],zoom_bounds[1][0]], zoom_bounds[0][1]-zoom_bounds[0][0], zoom_bounds[1][1]-zoom_bounds[1][0], fill=False,color='k',linewidth=2)
    pc = PatchCollection([rect],facecolor='None',edgecolor='k',linewidth=2)
    ax=plt.gca()
    ax.add_collection(pc)
    plt.savefig(os.path.join(savedir,saveprefix+'_black_zbox.png'))

    plt.figure()
    plt.imshow(tmp150,interpolation='none',cmap='bone')
    plt.clim(clim)
    plt.ylim(100,2800)
    plt.xticks([])
    plt.yticks([])
    if plot_coords:
        for ii,kk in enumerate(ra_tck_pos):
            plt.text(kk+50,2700,ra_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'))
        for ii,kk in enumerate(dec_tck_pos):
            plt.text(40,kk+200-ii*20,dec_tck_levels[ii],color='w',fontsize=12,bbox=dict(facecolor='black'),rotation=45-ii*6)
        plt.contour(ra,ra_contour_levels,colors='w', linestyles='dashed',linewidths=.75)
        plt.contour(dec,dec_contour_levels,colors='w',linestyles='dashed',linewidths=.75)
    plt.xlim(zoom_bounds[0])
    plt.ylim(zoom_bounds[1])
    print 'dec bounds'
    print dec[zoom_bounds[1][0],zoom_bounds[0][0]], dec[zoom_bounds[1][1],zoom_bounds[0][0]]
    print 'ra bounds'
    print ra[zoom_bounds[1][0],zoom_bounds[0][0]], ra[zoom_bounds[1][0],zoom_bounds[0][1]]
    plt.savefig(os.path.join(savedir,saveprefix+'_black_zboxzoom.png'))
 
    return ra, dec


