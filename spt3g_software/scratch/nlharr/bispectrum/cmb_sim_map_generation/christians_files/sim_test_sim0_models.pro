function beam_l,sigma_b,l
a=0.85
return,  (1.-a) * exp(-.5 * (.00292*l)^1.8) + a * exp(-.5 * (sigma_b*l)^1.5)
end
function beam_adj,dsigmas,fri,frj,l
 sigma_b0 = [0.00017172336, 0.00016751516]
 sigma_b = sigma_b0 + dsigmas
 
 return,  beam_l(sigma_b0(fri),l)*beam_l(sigma_b0(frj),l)/beam_l(sigma_b(fri),l)/beam_l(sigma_b(frj),l)
 
end

pro sim_test_sim0_models


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
iii = where(l eq 3000)

readcol,'/home/cr/paramfits/cosmomc.s10/camb/ml_20100826/ml_l10000_acc2_lensedCls.dat',l2,dl2,t1,t2,t3
dlcmb = [[dl2[l-2]],[dl2[l-2]],[dl2[l-2]]]

readcol,'/home/cr/paramfits/cosmomc.s10/camb/ml_20100915/ml_20100915_cmc_scalCls.dat',l2,dl2,t1,t2
ecl1 = dl2[2998]
ecl2 = dl2[3998]
dlcmb2=dlcmb*0
jj = where(l gt 4000,complement=cjj)
ljj=double(l(jj))
dlcmb2[cjj,0] = dl2[cjj+48]
dlcmb2[jj,0] = ecl2 * (ljj/4000d)^(alog(ecl2/ecl1) / alog(4000d/3000d))
dlcmb2[*,1]=dlcmb2[*,0]
dlcmb2[*,2]=dlcmb2[*,0]




llfac0 = l*(l+1d)/(2d*!dpi)
llfac = llfac0 / ((3000* 3001d)/(2d*!dpi))

czero_rg = 1.28
alpha_rg = -0.39
dlradio150 = czero_rg*llfac
dlradio220 = czero_rg * (219.3/151.7)^(2.*alpha_rg)*llfac
dlradio =[[dlradio150],[sqrt(dlradio150*dlradio220)],[dlradio220]]
print,czero_rg,1.28
print, czero_rg * (219.3/151.7)^(2.*alpha_rg), .95


czero_rg = 0.523872
dlradio150 = czero_rg*llfac
dlradio220 = czero_rg * (219.3/151.7)^(2.*alpha_rg)*llfac
dlradio2 =[[dlradio150],[sqrt(dlradio150*dlradio220)],[dlradio220]]



czero_dg = 9.5d
czero_cl= 4.5 
alpha_dg = 3.88
dg220 = (221.4/154.2)^(2.*alpha_dg)
dsfg150 = llfac*czero_dg
dldsfg = [[dsfg150],[dsfg150*sqrt(dg220)],[dsfg150*dg220 ]]
print,dldsfg[iii,*],(llfac0[iii])[0] *[6.6e-6,2.24439e-05,7.73299e-05]

readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/clustered_150.dat',l2,cclus
dclus = dlcmb*0.0
dlclus = cclus[48:*]*llfac0
dlclus /=dlclus[2950]
dlclus *=czero_cl
temp220 = (l/3000.)^(-.0000184825*2.*(221.4-154.2))
dlclus220 = dlclus * temp220 * (221.4/154.2)^(2*alpha_dg)
dclus[*,0]=dlclus
dclus[*,1]=sqrt(dlclus220*dlclus)
dclus[*,2]=dlclus220



czero_dg =9.24762 
czero_cl= 8.32451
alpha_dg = 4.18220
dg220 = (221.4/154.2)^(2.*alpha_dg)
dsfg150 = llfac*czero_dg
dldsfg2 = [[dsfg150],[dsfg150*sqrt(dg220)],[dsfg150*dg220 ]]
dclus2 = dlcmb*0.0
dlclus = cclus[48:*]*llfac0
dlclus /=dlclus[2950]
dlclus *=czero_cl
dlclus220 = dlclus * temp220 * (221.4/154.2)^(2*alpha_dg)
dclus2[*,0]=dlclus
dclus2[*,1]=sqrt(dlclus220*dlclus)
dclus2[*,2]=dlclus220





readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/tsz_sehgal.dat',l2,dtsz
dtszs = dlcmb*0.0
dtszs[*,0]=[fltarr(50),dtsz*llfac0]
dtszs[*,0]/=dtszs[iii,0] 
dtszs[*,0]*= 4.3
dtszs2 = dlcmb*0.0
dtszs2[*,0]=[fltarr(50),dtsz*llfac0]
dtszs2[*,0]/=dtszs[iii,0]
dtszs2[*,0]*= 1.80548



readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/ksz_sehgal.dat',l2,dksz
dkszs = dlcmb*0.0
dkszs[*,0]=[dksz[48:*]*llfac0]
dkszs[*,1]=[dksz[48:*]*llfac0]
dkszs[*,2]=[dksz[48:*]*llfac0]

dlth = dlradio + dkszs + dtszs + dlcmb + dldsfg +dclus

dlth2 = dlradio2 + dkszs + dtszs2 + dlcmb2 + dldsfg2 +dclus2


wins = fltarr(9951,45)
for i=1,45 do begin
    readcol,'~cr/s10/window_'+strtrim(string(i),2),l3,w3
    wins[*,i-1]=w3[0:9950]
endfor

dltmp = fltarr(9951,3)
dltmp2 = fltarr(9951,3)

a150 = 1 + randomn(seed) *.036
rX = 1 + randomn(seed)*.062

calfacs = a150^2 * [1.,rx,rx*rx]

dsigma150 = randomn(seed)*2.86e-6
dsigma220 = randomn(seed)*6.83e-6

dsigmas = [dsigma150,dsigma220]


dltmp[*,0] = dlth[*,0] * calfacs[0]* beam_adj(dsigmas,0,0,l)
dltmp[*,1] = dlth[*,1] * calfacs[1]* beam_adj(dsigmas,0,1,l)
dltmp[*,2] = dlth[*,2] * calfacs[2]* beam_adj(dsigmas,1,1,l)

bs = transpose(wins[*,0:14]) # dltmp
bs[*,1] = transpose(wins[*,15:29]) # dltmp[*,1]
bs[*,2] = transpose(wins[*,30:44]) # dltmp[*,2]
bs = reform(bs,45)

readcol,'~cr/s10/covariance_150_cross_220_s10.txt',cov
cov = reform(cov,45,45)
evec=cov
TRIRED, evec, D, E
; Compute the eigenvalues (returned in vector D) and                                       
; the eigenvectors (returned in the rows of the array A):                                  
TRIQL, D, E, evec
err = sqrt(d)

berrdiag =  randomn(seed,45)*err
berr = evec * berrdiag
bs0=bs+berr


ccal1 = 0.998411
ccal2 = 1.04879
dsig1 = 1.61746e-06
dsig2 =  -1.32045e-05
dsigs2 = [dsig1,dsig2]
dltmp2[*,0] = dlth2[*,0] * ccal1*ccal1* beam_adj(dsigs2,0,0,l)
dltmp2[*,1] = dlth2[*,1] * ccal1*ccal2* beam_adj(dsigs2,0,1,l)
dltmp2[*,2] = dlth2[*,2] * ccal2*ccal2* beam_adj(dsigs2,1,1,l)
bs2 = transpose(wins[*,0:14]) # dltmp2
bs2 = reform(bs2,45)

readcol,'~cr/s10_mf_test/v1/test_0',bb,bsin

plot,bs
oplot,bs2,linestyle=2
oplot,bsin,color=200
oplot,bs0,color=200,linestyle=2
stop
end
