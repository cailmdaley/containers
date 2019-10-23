### A tool for theoretical noise estimation
#### Component values are based on spec sheets and the current ICE boards
#### Formulas come from the noise memo, there may be some minor issues, but are likely pretty close

k_b=1.38e-23  #boltzmann
T=300  #K
import numpy as np
import matplotlib.pyplot as plt
import pdb

def mezz_transfer_function(rf1=300.,rf2=200.,rg=300.):
    '''
    RF1: feedback resistor on the first stage transimpedance amplifier
    RF2: feedback resistor on the second stage fully differential amplifier
    Rg : variable gain resistor before second stage amplifier
    gain= mezz_transfer_function(rf1=300,rf2=200,rg=300,r3 =50)
    voltage_out= gain*current_from_DAC
    '''
    gain= rf1*rf2/rg
    return gain
    

def carrier_transfer_function(rf1=300.,rf2=200.,rg=300.,r3=50.,r4=20.,rwireharness=0.,rbias=0.03,rtes=1.5):
    '''
    R3: resistors on mezzanine after amplifiers
    r4: resistors on squid controller board
    rwireharness: resistance in series due to wire harness
    rbias: bias resistor
    rtes: tes resistance

    NB: this assumes that the LC circuit has zero impedance 
    '''
    
    z_total=(rbias*rtes)/(rbias+rtes)
    mezz_gain=mezz_transfer_function(rf1=rf1,rf2=rf2,rg=rg)
    carrier_tf=mezz_gain/(2*r3+4*r4+rwireharness+z_total)*z_total    
    return carrier_tf

def nuller_transfer_function(rf1=300.,rf2=200.,rg=300.,r3=50., r5=100., r6=750.):
    '''
    r3: resistors on mezzanine after amplifiers
    r5: resistor across the output lines at input of squid controller
    56: resistors on squid controller board
    lsq: squid input coil inductance
    '''
    
    mezz_gain=mezz_transfer_function(rf1=rf1,rf2=rf2,rg=rg)
    zthev=(r5*4*r6)/(r5+4.*r6)
    nuller_tf=mezz_gain/(2*r3+zthev)*r5/(4.*r6+r5)
    return nuller_tf


def calc_carrier_noise(rf1=300., rf2=200., rg=300., r2=250., r3=50., r4=20., rwireharness=0., rbias=0.03, rtes=1.5):
    # DAC noise
    carrier_tf=carrier_transfer_function(rf1=rf1,rg=rg, rf2=rf2,r3=r3,r4=r4,rwireharness=rwireharness,rbias=rbias,rtes=rtes)
    e_dac=50e-12*carrier_tf
    #print(e_dac)
    #first stage amp
    z_total=(rbias*rtes)/(rbias+rtes)
    e_lt6321=calc_lt6321_noise(rf1=rf1)*rf2/rg/(2*r3+4*r4+rwireharness+z_total)*z_total
    #print(e_lt6321)
    # second stage amp
    e_ths4131=calc_ths4131_noise(rf2=rf2,rg=rg,r2=r2) /(2*r3+4*r4+rwireharness+z_total)*z_total
    #print(e_ths4131)
    # mezz resistors
    e_mezz_res=np.sqrt(4*k_b*T*(2*r3)) /(2*r3+4*r4+rwireharness+z_total)*z_total
    #print(e_mezz_res)
    #squid controller board resistors ( & wire harness, which is not quite right because T not 300)
    e_sqbd_res=np.sqrt(4*k_b*T*(4*r4+rwireharness)) /(2*r3+4*r4+rwireharness+z_total)*z_total
    #print(e_sqbd_res)
    e_total= np.sqrt( e_dac**2+e_lt6321**2 + e_ths4131**2+e_mezz_res**2+e_sqbd_res**2)
    i_total = e_total/rtes
    return e_total, i_total

def calc_nuller_noise(rf1=300., rf2=200., rg=300., r2=250., r3=50.,r5=100., r6=750.):

    # dac noise
    nuller_tf=nuller_transfer_function(rf1=rf1,rf2=rf2,rg=rg,r3=r3, r5=r5, r6=r6)
    i_dac=50e-12*nuller_tf
    #print(i_dac)
    #first stage amp
    zthev=(r5*4*r6)/(r5+4.*r6)
    i_lt6321=calc_lt6321_noise(rf1=rf1)*rf2/rg/(2*r3+zthev)*r5/(4.*r6+r5)
    #print(i_lt6321)
    #second stage amp
    i_ths4131=calc_ths4131_noise(rf2=rf2,rg=rg,r2=r2)/(2*r3+zthev)*r5/(4.*r6+r5)
    #print(i_ths4131)
    #mezz res
    i_mezz_res=np.sqrt(4*k_b*T*(2*r3)) /(2*r3+zthev)*r5/(4.*r6+r5)
    #print(i_mezz_res)
    #sqbd_res
    i_sqbd_res1=np.sqrt(4*k_b*T*(r5))/(4.*r6+r5)
    #print(i_sqbd_res1)
    #sqbd res
    i_sqbd_res2=np.sqrt(4*k_b*T*(4*r6))/(4.*r6+r5)
    #print(i_sqbd_res2)
    
    i_total=np.sqrt(i_dac**2+i_lt6321**2+i_ths4131**2+i_mezz_res**2+i_sqbd_res1**2+i_sqbd_res2**2)    
    return i_total

def calc_demod_noise(r1d=10.,rf1d=150., r2d=100., rf2d=400., r3d=50.,r4d=42., r5d=100., rdynsq=500.,zsquid=500.):
    #first stage amplifier
    e_lt6200=calc_lt6200_noise(rf1d=rf1d,r1d=r1d,rdynsq=rdynsq)/((rf1d+r1d)/r1d)
    #print(e_lt6200)
    # second stage amplifier
    e_ths4131_1=calc_ths4131_noise(rf2=rf2d, rg=r2d, r2=0)/((rf1d+r1d)/r1d)/(rf2d/r2d)
    #print(e_ths4131_1)
    # third stage amplifier
    e_ths4131_2=calc_ths4131_noise(rf2=rf2d, rg=r2d, r2=0)/((rf1d+r1d)/r1d)/(rf2d/r2d)/(rf2d/r2d)
    #print(e_ths4131_2)
    e_total=np.sqrt(e_lt6200**2+e_ths4131_1**2+e_ths4131_2**2)
    #print(e_total)
    i_total=e_total/zsquid
    return i_total

def calc_lt6321_noise(rf1=300.):
    e_n_intrinsic=1.1e-9
    i_n_intrinsic=2.2e-12
    k_b=1.38e-23  #boltzmann
    T=300  #K
    e_total= np.sqrt(e_n_intrinsic**2+(i_n_intrinsic*rf1)**2+4*k_b*T*rf1)
    e_total=e_total*np.sqrt(2)  # factor of sqrt of two because the calc is for a single amp, but the configuration of this is actually a double (identical features)
    return e_total

def calc_ths4131_noise(rf2=200., rg=300., r2=250.):
    e_n_intrinsic=1.3e-9
    i_n_intrinsic=1e-12
    k_b=1.38e-23  #boltzmann
    T=300  #K
    e_johnson=np.sqrt(2)*np.sqrt(  4*k_b*T*rf2 + 4*k_b*T*rg*(rf2/rg)**2 +4*k_b*T*r2*((rf2+rg)/rg)**2)
    e_internal=np.sqrt(  ( (rf2+rg)/rg*e_n_intrinsic)**2 +  2* ( i_n_intrinsic* (rf2+r2+ rf2*r2/rg))**2)
    e_total=np.sqrt(e_johnson**2+e_internal**2)
    return e_total

def calc_lt6200_noise(r1d=10., rf1d=150., rdynsq=500.):
    e_n_intrinsic=1.1e-9
    i_n_intrinsic=2.2e-12
    e_johnson= np.sqrt(4*k_b*T*rf1d + 4*k_b*T*r1d*((rf1d+r1d)/r1d)**2)
    e_internal= np.sqrt( e_n_intrinsic**2*((rf1d+r1d)/r1d)**2 +(i_n_intrinsic*rdynsq*(rf1d+r1d)/r1d)**2+(i_n_intrinsic*rf1d)**2)
    e_total=np.sqrt(e_johnson**2+e_internal**2)
    return e_total

def calc_bias_resistor_noise(rbias=0.03,rtes=1.5,Tbias=4):
    #T=4  #cold bias resistor
    e_bias=np.sqrt(4*k_b*Tbias*rbias)
    i_bias=e_bias/rtes
    return i_bias

def calc_squid_noise(squid_noise=4.5e-12):
    return squid_noise

def calc_current_sharing(freq=1.e6,L_squid=60.e-9,rtes=1.5,Lstrip=46e-9):
    amp_fac=np.abs(1.+ 2j*np.pi*freq*L_squid/(rtes+2j*np.pi*freq*Lstrip))
    return amp_fac

def calc_total_noise_current(rf1=300., rf2=200., rg=300., r2=250., r3=50., r4=20., rwireharness=0., rbias=0.03, rtes=1.5,r5=100., r6=750.,r1d=10.,rf1d=150., r2d=100., rf2d=400., r3d=50.,r4d=42., r5d=100., rdynsq=500.,zsquid=500.,L_squid=60e-9,Lstrip=46e-9,current_share_factor=False,rtes_v_f=False,freqs=[0.],Tbias=4.,add_excess=[],add_excess_demod=[]):
    #note this is the noise current through the squid coil
    e_total_c,i_total_c=calc_carrier_noise(rf1=rf1, rf2=rf2, rg=rg, r2=r2, r3=r3, r4=r4, rwireharness=0., rbias=rbias, rtes=rtes)
    #print(i_total_c*1e12)
    i_total_n=calc_nuller_noise(rf1=rf1, rf2=rf2, rg=rg, r2=r2, r3=r3,r5=r5, r6=r6)
    #print(i_total_n*1e12)
    i_total_d=calc_demod_noise(r1d=r1d,rf1d=rf1d, r2d=r2d, rf2d=rf2d, r3d=r3d,r4d=r4d, r5d=r5d, rdynsq=rdynsq,zsquid=zsquid)
    #print(i_total_d*1e12)
    i_bias=calc_bias_resistor_noise(rbias=rbias,rtes=rtes,Tbias=Tbias)
    #print(i_bias*1e12)
    i_squid=calc_squid_noise()
    #print(i_squid*1e12)
    i_total=np.sqrt(i_total_c**2+i_total_n**2+i_total_d**2+i_bias**2+i_squid**2)
    if len(freqs)==1 and freqs[0]==0.:
        freqs=np.arange(1e6,6e6,1e4)
        i_total=i_total*(np.zeros_like(freqs)+1)
    if current_share_factor:
        if rtes_v_f:
            #hardcoded 20% decrease in rtes as a function of frequency
            rtes_ch=rtes*np.arange(1,0.6,-0.4/len(freqs))
            amp_fac=np.zeros(len(freqs))
            for ii,kk in enumerate(freqs):
                amp_fac[ii]=calc_current_sharing(freq=kk,L_squid=L_squid,Lstrip=Lstrip,rtes=rtes_ch[ii])            
        else:
            amp_fac=calc_current_sharing(freq=freqs,L_squid=L_squid,Lstrip=Lstrip,rtes=rtes)
            i_total_d_cs=amp_fac*i_total_d
        i_total=np.sqrt(i_total_c**2+i_total_n**2+i_total_d_cs**2+i_bias**2+i_squid**2)
    if len(add_excess_demod)>0:
        i_total=np.sqrt(i_total_c**2+i_total_n**2+(i_total_d_cs*add_excess_demod)**2+i_bias**2+i_squid**2)
    if len(add_excess)>0:
        i_total=i_total*add_excess
    
    return freqs,i_total

    

def calc_total_noise_power(rf1=300., rf2=200., rg=300., r2=250., r3=50., r4=20., rwireharness=0., rbias=0.03, rtes=1.5,r5=100., r6=750.,r1d=10.,rf1d=150., r2d=100., rf2d=400., r3d=50.,r4d=42., r5d=100., rdynsq=500.,zsquid=500.,L_squid=60e-9,Lstrip=46e-9,pelectrical=5e-12,freqs=[0.], current_share_factor=False,Tbias=4.,add_excess=[], rtes_v_f=False,add_excess_demod=[]):
    '''
    rtes: resistance of the TES
    rg: gain resistor on mezzanine that controls analog gain (300 = gain 0, 100 = gain 15)
    z_squid: squid gtransimpedance
    rdynsq: squid dynamic  impedance
    pelectrical: electrical bias power assumed, used for calculating vbias
    L_squid: input inductance of the SQUID
    Lstrip: inductance of the stripline (or parasitic between comb summing point and SQUID input coil
    freqs: bias frequency to perform calculation of current sharing at.  If 0., will be 10kHz points between 1 & 6 MHz
    current_share_factor: True = turn on current sharing
    add_excess: <number> = add arbitrary excess factor (False = off)
    Tbias: operating temperature of the bias resistor
    '''

    
    freqs,i_total=calc_total_noise_current(rf1=rf1, rf2=rf2, rg=rg, r2=r2, r3=r3, r4=r4, rwireharness=0., rbias=rbias, rtes=rtes,r5=r5, r6=r6,r1d=r1d,rf1d=rf1d, r2d=r2d, rf2d=rf2d, r3d=r3d,r4d=r4d, r5d=r5d, rdynsq=rdynsq,zsquid=zsquid,L_squid=L_squid,Lstrip=Lstrip,current_share_factor=current_share_factor, rtes_v_f=rtes_v_f,freqs=freqs,Tbias=Tbias,add_excess=add_excess,add_excess_demod=add_excess_demod)

    vbias=np.sqrt(pelectrical*rtes)
    #print('Vbias: %5.2f \mu Volts' % (vbias*1e6))
    #note the extra sqrt(2) for AC biased system
    nep_readout=vbias*i_total/np.sqrt(2)

    return freqs,nep_readout
    

def plot_nep_vs_rtes(rtes_range=np.arange(0.25,2.25,0.25), zsquid=500., rg=300, pelectrical=3e-12,L_squid=60e-9, Lstrip=46e-9,current_share_factor=False):
    nep_arr=np.zeros_like(rtes_range)
    
    plt.figure()
    plp=[]
    legstr=[]
    for ii,kk in enumerate(rtes_range):
        freqs,nep_arr=calc_total_noise_power(rtes=kk,rg=rg, zsquid=zsquid,pelectrical=pelectrical,current_share_factor=current_share_factor,L_squid=L_squid,Lstrip=Lstrip)
        p1,=plt.plot(freqs/1e6,nep_arr*1e18)
        plp.append(p1)
        legstr.append(str(kk)+'$\Omega$')
    plt.xlabel('Bias Frequency (MHz)')
    plt.ylabel('Readout Noise Equivalent Power (aW/$\sqrt{Hz}$)')
    plt.legend(plp,legstr)
    plt.title('$P_{electrical}$ = %5.2f pW' % (pelectrical*1e12))
    return

def plot_2dgrid_nep_vs_f_rtes(rtes_range=np.arange(0.25,2.25,0.0025), zsquid=500., rg=300, pelectrical=3e-12,L_squid=60e-9, Lstrip=46e-9,current_share_factor=False,freqs=[0.],add_excess=[],Tbias=4.,label='',add_excess_demod=[]):        
    #freqs=np.arange(1.e6,6.0e6,1.e4)
    #print(len(freqs))
    #print(len(rtes_range))
    nep_2darr=np.zeros((len(rtes_range),len(freqs)))
    for ii,kk in enumerate(rtes_range):
        freqs_tmp,nep_tmp=calc_total_noise_power(rtes=kk,rg=rg, zsquid=zsquid,pelectrical=pelectrical,current_share_factor=current_share_factor,L_squid=L_squid,freqs=freqs,Lstrip=Lstrip,add_excess=add_excess,Tbias=Tbias,add_excess_demod=[])
        nep_2darr[ii,:]=nep_tmp
        
    plt.figure()
    plt.imshow(nep_2darr*1e18,interpolation='none',origin='lower')
    plt.xlabel('Bias Frequency (MHz)')
    plt.ylabel('TES Resistance ($\Omega$)')
    plt.colorbar(label='NEP Readout (aW/$\sqrt{Hz}$)')
    if label:
        plt.title(label)
    else:
        plt.title('$P_{electrical}$ = %5.2f pW' % (pelectrical*1e12))
    plt.xticks(np.arange(0,len(freqs),100), freqs[np.arange(0,len(freqs),100)]/1e6)
    plt.yticks(np.arange(0,len(rtes_range),np.int32(np.round(0.25/np.diff(rtes_range)[0]))), np.round(rtes_range[np.arange(0,len(rtes_range),np.int32(np.round(0.25/np.diff(rtes_range)[0])))],2))
    return 

def plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=np.arange(0.25,2.25,0.0025), zsquid=500., rg=300, photon_noise=0, phonon_noise=0, pelectrical=3e-12,L_squid=60e-9, Lstrip=46e-9,current_share_factor=False,freqs=[0.],add_excess=[],Tbias=4,label='',add_excess_demod=[]):
    #freqs=np.arange(1.e6,6.01e6,1.e4)
    #print(len(freqs))
    #print(len(rtes_range))
    nep_2darr=np.zeros((len(rtes_range),len(freqs)))
    for ii,kk in enumerate(rtes_range):
        freqs_tmp,nep_tmp=calc_total_noise_power(rtes=kk,rg=rg, zsquid=zsquid,pelectrical=pelectrical,current_share_factor=current_share_factor,L_squid=L_squid,freqs=freqs,Lstrip=Lstrip,add_excess=add_excess,Tbias=Tbias,add_excess_demod=[])
        nep_2darr[ii,:]=nep_tmp
    nep_det_only=np.sqrt(phonon_noise**2 + photon_noise**2)
    nep_total_arr=np.sqrt(nep_2darr**2+phonon_noise**2 + photon_noise**2)
    nep_ratio_arr=nep_total_arr/nep_det_only
    plt.figure()
    plt.imshow(nep_ratio_arr,interpolation='none',origin='lower')
    plt.xlabel('Bias Frequency (MHz)')
    plt.ylabel('TES Resistance ($\Omega$)')
    #    if len(add_excess)>0:
    #        plt.clim(1.0,1.8)
    #    else:
    #        plt.clim(1.0,1.2)
    plt.colorbar(label='Increase in Noise from Readout')
    ctrlevels=np.arange(1.05,1.4,0.05)
    ctr=plt.contour(np.arange(len(freqs)),np.arange(len(rtes_range)),nep_ratio_arr,levels=ctrlevels,colors='w')
    plt.clabel(ctr,inline=1,fontsize=10)
    if label:
        plt.title(label)
    else:
        plt.title('$P_{electrical}$ = %5.2f pW' % (pelectrical*1e12))
    plt.xticks(np.arange(0,len(freqs),100), freqs[np.arange(0,len(freqs),100)]/1e6)
    plt.yticks(np.arange(0,len(rtes_range),np.int32(np.round(0.25/np.diff(rtes_range)[0]))), np.round(rtes_range[np.arange(0,len(rtes_range),np.int32(np.round(0.25/np.diff(rtes_range)[0])))],2))
    return 

    

def plot_s4_nep_grids(spt3g_theory=False,spt3g_scaled_theory=False, spt3g_lowII_theory=False, spt3g_lowII_scaled_theory=False,mkSQUID_theory=False,mkSQUID_scaled_theory=False,mkSQUID_lowII_theory=False,mkSQUID_lowII_scaled_theory=False,all=False):
    if all:
        spt3g_theory=True
        spt3g_scaled_theory=True
        spt3g_lowII_theory=True
        spt3g_lowII_scaled_theory=True
        mkSQUID_theory=True
        mkSQUID_scaled_theory=True
        mkSQUID_lowII_theory=True
        mkSQUID_lowII_scaled_theory=True
    
    parameter_dict={}
    parameter_dict['lat']={}
    parameter_dict['sat']={}

    #LAT parameters
    parameter_dict['lat']['pelectricals']=np.array([0.5e-12,0.4e-12,2.2e-12,2.4e-12,5.0e-12,17.2e-12,20.2e-12])
    parameter_dict['lat']['photon_noise']=np.array([5.17,4,14.8,15.1,28.3,71.7,91])*1e-18
    parameter_dict['lat']['phonon_noise']=np.array([3.86,3.3,8.4,8.6,12.5,23.2,25.1])*1e-18
    parameter_dict['lat']['labels']=['lat20ghz','lat27ghz','lat39ghz','lat93ghz','lat145ghz','lat225ghz','lat278ghz']
    
    #SAT parameters
    parameter_dict['sat']['pelectricals']=np.array([0.8,1.7,3.3,4.7,3.0,5.0,9.5,13.1])*1e-12
    parameter_dict['sat']['photon_noise']=np.array([8.5,15,25,31,23,32,57,72])*1e-18
    parameter_dict['sat']['phonon_noise']=np.array([4.6,6.5,10,13,10,14,20,23])*1e-18
    parameter_dict['sat']['labels']=['sat30ghz','sat40ghz','sat85ghz','sat145ghz','sat95ghz','sat155ghz','sat220ghz','sat270ghz']    
    
    #SPT-3G excess scaling
    freqs=np.arange(1e6,6.0e6,1e4)
    ff=np.arange(0,len(freqs))
    excess_fac=np.polyval([5e-6,0,1.45],ff)

    #basic inputs
    rtes_range=np.arange(0.25,2.25,0.005)
    zsquid=500.
    rg=300.
    L_squid=60e-9
    Lstrip=46e-9
    current_share_factor=True
    
    if spt3g_theory:
        lbltxt='SPT3G theory'
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc1_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc1_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc1_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc1_'+parameter_dict[jj]['labels'][ii]+'.png')

    if spt3g_scaled_theory:
        lbltxt='SPT3G scaled'
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs,add_excess=excess_fac,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc2_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs,add_excess=excess_fac,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc2_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc1_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc1_'+parameter_dict[jj]['labels'][ii]+'.png')

    if spt3g_lowII_theory:
        lbltxt='Lsquid=10nH theory'
        #Lower input inductance squid coil.  No excess scaling
        L_squid=10e-9
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc3_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc3_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc3_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc3_'+parameter_dict[jj]['labels'][ii]+'.png')
        L_squid=60e-9  #reset the value
    if spt3g_lowII_scaled_theory:
        lbltxt = 'Lsquid=10nH scaled'
        #Lower input inductance squid coil.  No excess scaling
        L_squid=10e-9
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs, add_excess=excess_fac,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc4_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs, add_excess=excess_fac,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc4_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc4_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc4_'+parameter_dict[jj]['labels'][ii]+'.png')
        L_squid=60e-9  #reset the value

    if mkSQUID_theory:
        lbltxt='mK SQUID, Lsquid=60nH theory'
        Lstrip=8e-9
        Tbias=0.25
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc5_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc5_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc5_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc5_'+parameter_dict[jj]['labels'][ii]+'.png')
        Lstrip=46e-9
        Tbias=4
        
    if mkSQUID_scaled_theory:
        lbltxt='mK SQUID, Lsquid=60nH scaled'
        Lstrip=8e-9
        Tbias=0.25
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs,add_excess=excess_fac,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc6_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs,add_excess=excess_fac,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc6_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc6_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc6_'+parameter_dict[jj]['labels'][ii]+'.png')
        Lstrip=46e-9
        Tbias=4

    if mkSQUID_lowII_theory:
        lbltxt='mK SQUID, Lsquid=10nH theory'
        Lstrip=8e-9
        L_squid=10e-9
        Tbias=0.25
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc7_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc7_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc7_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc7_'+parameter_dict[jj]['labels'][ii]+'.png')
        Lstrip=46e-9
        L_squid=60e-9
        Tbias=4
        
    if mkSQUID_lowII_scaled_theory:
        lbltxt='mK SQUID, Lsquid=10nH scaled'
        Lstrip=8e-9
        L_squid=10e-9
        Tbias=0.25
        for jj in parameter_dict.keys():
            for ii,kk in enumerate(parameter_dict[jj]['pelectricals']):
                plot_2dgrid_nep_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,freqs=freqs,add_excess=excess_fac,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_vs_f_rtes_sc8_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
                plot_2dgrid_totalnep_increase_vs_f_rtes(rtes_range=rtes_range, zsquid=zsquid, rg=rg, pelectrical=kk,L_squid=L_squid, Lstrip=Lstrip,current_share_factor=current_share_factor,photon_noise=parameter_dict[jj]['photon_noise'][ii],phonon_noise=parameter_dict[jj]['phonon_noise'][ii],freqs=freqs,add_excess=excess_fac,Tbias=Tbias,label=lbltxt)
                plt.tight_layout()
                plt.savefig('/Users/benderan/Desktop/lowR_dfmux/noise_memo/s4memo_figures/'+'nepreadout_increase_vs_f_rtes_sc8_'+parameter_dict[jj]['labels'][ii]+'.png')
                plt.close()
            for ii in range(len(parameter_dict[jj]['pelectricals'])):
                print('nepreadout_vs_f_rtes_sc8_'+parameter_dict[jj]['labels'][ii]+'.png')
                print('nepreadout_increase_vs_f_rtes_sc8_'+parameter_dict[jj]['labels'][ii]+'.png')
        Lstrip=46e-9
        Tbias=4
        L_squid=60e-9

    return

def plot_total_nep_vs_rtes(rtes_range=np.arange(0.25,2.25,0.25), zsquid=500., rg=300, pelectrical=3e-12,current_share_factor=False):
    nep_arr=np.zeros_like(rtes_range)
    plp=[]
    photon_nep=[33e-18, 39e-18, 46e-18, 53e-18]
    photon_nep=[]
    for kk in pelectrical:
        photon_nep.append(calc_photon_noise(poptical=kk,dnu=30e9,nu0=150e9,eps=1.))        
    plt.figure()
    for jj,pp in enumerate(pelectrical):
        for ii,kk in enumerate(rtes_range):
            nep_arr[ii]=calc_total_noise_power(rtes=kk,rg=rg, zsquid=zsquid,pelectrical=pp,excess_factor=excess_factor)
        p1,=plt.plot(rtes_range,np.sqrt(nep_arr**2+(photon_nep[jj])**2)/photon_nep[jj],'.',label='P$_{opt}$ = P$_{elec}$ = '+np.str(np.round(pp*1e12))+'pW')
        plp.append(p1)
    plt.xlabel('TES Resistance ($\Omega$)')
    plt.ylabel('Fractional Increase in Total NEP from Readout')
    
    
    plt.legend()
    return
    




def calc_photon_noise(poptical=7.74e-12, dnu=30e9,nu0=150e9,eps=1.):
    h=6.62607004e-34
    popt=np.sqrt(2*h*nu0*poptical)
    print(popt*1e18)
    #for a polarization sensitive detector
    pcorr=np.sqrt(2.*eps*poptical**2/dnu)
    print(pcorr*1e18)
    ptotal=np.sqrt(popt**2+pcorr**2)
    print(ptotal*1e18)
    return ptotal

#NB. A quick check suggests that Matt & I agree on the photon noise calculation

def calc_phonon_noise():
    print('test')
    return

def calc_johnson_noise():
    return


#Dobbs noise sim gets responsivity from the sim and uses that instead of vbias to translate to NEP
#  1/corrected_vbias*np.sqrt(2)*(L/L+1)  = responsivity = dI/dP
# NEP = NEI/ (dI/dP)
