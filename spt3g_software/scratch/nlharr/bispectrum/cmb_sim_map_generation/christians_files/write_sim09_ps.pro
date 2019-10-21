pro write_sim09_ps,freq

pi=!dpi



calc=dblarr(25000)
calccmb=dblarr(25000)
ll=dindgen(25000)+2
dlfac = ll*(ll+1.)/(2*!pi)
fwhm2sig = 1./sqrt(8.*alog(2.))
   

calc = dlfac / ((3000.*3001.)/(2.*!pi))
                                ; I've normalized with respect to
                                ; D3000=1
calcell = calc /( ll/3000.) ; likewise D3000=1


facs = estimate_pointsources_sim09()
matrix = [[facs.ps[0],facs.ps[3],facs.ps[4]],$
          [facs.ps[3],facs.ps[1],facs.ps[5]],$
          [facs.ps[4],facs.ps[5],facs.ps[2]]]
ellmatrix = [[facs.ellclus[0],facs.ellclus[3],facs.ellclus[4]],$
             [facs.ellclus[3],facs.ellclus[1],facs.ellclus[5]],$
             [facs.ellclus[4],facs.ellclus[5],facs.ellclus[2]]]
evecs=matrix
trired,evecs,evals,e
triql,evals,e,evecs

ii0=1
if fix(freq) eq 220 then ii0=2
if fix(freq) eq 90 then ii0=0
amps = sqrt(evals)
effs=fltarr(3)
effs[0] = amps[0]*evecs[ii0,0]
effs[1] = amps[1]*evecs[ii0,1]
effs[2] = amps[2]*evecs[ii0,2]

aps = (effs[1]^2 + effs[2]^2)
ps = calc*aps

;check if planned N makes sense
jj=abs(evals)
ind = where(jj gt 1e-4*max(jj),nbig,complement=cnd)
if nbig ne 2 then stop
if nbig lt 3 then effs[cnd]=0

evecs=ellmatrix
trired,evecs,evals,e
triql,evals,e,evecs

jj=abs(evals)
ind = where(jj gt 1e-4*max(jj),nbig,complement=cnd)
if nbig ne 1 then stop
ampsclus = sqrt(evals)

effsclus=fltarr(3)
effsclus[0] = ampsclus[0]*evecs[ii0,0]
effsclus[1] = ampsclus[1]*evecs[ii0,1]
effsclus[2] = ampsclus[2]*evecs[ii0,2]
if nbig lt 3 then effsclus[cnd]=0
;stop
aclus = effsclus[2]^2
psclus = aclus * calcell
sfreq = strtrim(string(freq),2)

a = read_ascii('~cr/cmb_models/foreground_sim09_'+sfreq+'.txt')
psold = reform(a.field1[1,*])
stop
end




