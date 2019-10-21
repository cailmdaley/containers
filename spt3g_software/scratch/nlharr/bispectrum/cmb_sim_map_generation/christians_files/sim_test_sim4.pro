function beam_l,sigma_b,l
a=0.85
return,  (1.-a) * exp(-.5 * (.00292*l)^1.8) + a * exp(-.5 * (sigma_b*l)^1.5)
end
function beam_adj,dsigmas,fri,frj,l
 sigma_b0 = [0.00017172336, 0.00016751516]
 sigma_b = sigma_b0 + dsigmas
 
 return,  beam_l(sigma_b0(fri),l)*beam_l(sigma_b0(frj),l)/beam_l(sigma_b(fri),l)/beam_l(sigma_b(frj),l)
 
end

pro sim_test_sim4

nsim = 10

seed=7
;

;    ! Effective frequencies for the bands                                                             
;   ! The order is (dutsy clustered,
                                ;   dusty poisson, radio poisson, ksz,
                                ;   tsz)                          
;  real, dimension(5,2), parameter :: spt_eff_fr \
   ;     = (/ (/ 154.2, 154.2, 151.7, 150, 153/), \
 ;            (/ 221.4, 221.4, 219.3, 220, 220.2/) /)


l = dindgen(9951)+50


readcol,'/home/cr/paramfits/cosmomc.s10/camb/ml_20100826/ml_l10000_acc2_lensedCls.dat',l2,dl2,t1,t2,t3
dlcmb = [[dl2[l-2]],[dl2[l-2]],[dl2[l-2]]]

llfac0 = l*(l+1d)/(2d*!dpi)
llfac = llfac0 / ((3000* 3001d)/(2d*!dpi))
dlradio150 = 1.28*llfac
dlradio220 = 0.95*llfac
dlradio =[[dlradio150],[sqrt(dlradio150*dlradio220)],[dlradio220]]

dldsfg = [[llfac0*6.6e-6],[llfac0* 2.24439e-05],[llfac0*7.73299e-05 ]]

readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/tsz_sehgal.dat',l2,dtsz
dtszs = dlcmb*0.0
dtszs[*,0]=[fltarr(50),dtsz*llfac0]*0.5

readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/ksz_sehgal.dat',l2,dksz
dkszs = dlcmb*0.0
dkszs[*,0]=[dksz[48:*]*llfac0]
dkszs[*,1]=[dksz[48:*]*llfac0]
dkszs[*,2]=[dksz[48:*]*llfac0]

dlth = dlradio + dkszs + dtszs + dlcmb + dldsfg


readcol,'/home/cr/paramfits/cosmomc.s10/sim4/suxp_lma10000_spt4',ls,dl1,dl2,dl3,clcmb0,clf1,clf2,clf3

dlsfac = ls*(ls+1d)/(2d*!dpi)
dlcmb0=clcmb0*dlsfac
dlf1=clf1*dlsfac
dlf2=clf2*dlsfac
dlf3=clf3*dlsfac

inds = where(l ge 1500)
dlfm =  dlradio + dkszs + dtszs + dldsfg
dlfm=dlfm[inds,*]
dlcmbm = dlcmb[inds,*]
 plot,ls,dlcmbm[*,0]-dlcmb0,xr=[2500,4500]
oplot,ls,dlcmbm[*,0]/1.036-dlcmb0,color=200


stop
end
