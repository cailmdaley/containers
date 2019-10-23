;;;
; NAME: get_highel_theory_dls.pro
; PURPOSE:
;   Get the input spectrum for the highel2500 sims in dl.
;
; INPUTS:
;   None
;
; NOTES:
; 1) expects file at scripts/tf_prep[run].sav
;
; MODIFICATION HISTORY:
;  07/31/2013: (KTS) Created
;;;

FUNCTION get_highel2500_theory_dls, cl=cl, remake=remake, fname=fname, stopit=stopit
    theoryspecdir = '/home/lizinvt/highel2500/cmb_models/'
    theoryspectrumfiles=[[theoryspecdir+'ml_20120911_lensedCls.dat',$
                          theoryspecdir+'dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake25000.txt', $
                          theoryspecdir+'dl_shaw_tsz_s10_97p6ghz_fake25000.txt',$
                          theoryspecdir+'poisson_90ghz.txt',$
                          theoryspecdir+'psclustered_90ghz.txt'],$
                         [theoryspecdir+'ml_20120911_lensedCls.dat',$
                          theoryspecdir+'dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake25000.txt',$
                          theoryspecdir+'dl_shaw_tsz_s10_152p9ghz_fake25000.txt',$
                          theoryspecdir+'poisson_150ghz.txt',$
                          theoryspecdir+'psclustered_150ghz.txt'],$
                         [theoryspecdir+'ml_20120911_lensedCls.dat',$
                          theoryspecdir+'dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake25000.txt',$
                          theoryspecdir+'dl_shaw_tsz_s10_218p1ghz_fake25000.txt',$
                          theoryspecdir+'poisson_220ghz.txt',$
                          theoryspecdir+'psclustered_220ghz.txt']]
    
    nfreq=3
    nell = 25000
    dl_th = dblarr(nfreq+1, nell)
    
    for ifreq=0, nfreq-1 do begin
        ; read in .dat file
        readcol, theoryspectrumfiles[0,ifreq], aell, att, aee, abb, ate
        aell = [0.0, 1.0, aell]
        att = [0.0, 0.0, att]
        ll=dindgen(nell)+2
        newatt = interpol([att,0.0], [aell,25001], [0.0, 1.0, ll])
        calccmb = ll*0.0
        calccmb[0:24999] = newatt[0:24999]
        
        ; ksz
        readcol, theoryspectrumfiles[1,ifreq], l1, y1
        dl1 = interpol(y1, l1, ll)
        
        ; tzs
        readcol, theoryspectrumfiles[2,ifreq], l2, y2
        dl2 = interpol(y2, l2, ll)
        
        ; poisson
        readcol, theoryspectrumfiles[3,ifreq], l3, y3
        dl3 = interpol(y3, l3, ll)
        
        ; clustered
        readcol, theoryspectrumfiles[4,ifreq], l4, y4
        dl4 = interpol(y4, l4, ll)
        
        ; combine
        dl_th[ifreq+1,*] = calccmb + dl1 + dl2 + dl3 + dl4
    endfor    
    dl_th[0,*] = ll

    spectrum = dl_th

    ; if requested, convert to Cl
    if keyword_set(cl) then begin
        cl_th = dl_th*0.
        cl_th[0,*] = dl_th[0,*]
        cl_th[1,*] = dl_th[1,*]*2*!DPI/(dl_th[0,*]*(dl_th[0,*]+1))
        cl_th[2,*] = dl_th[2,*]*2*!DPI/(dl_th[0,*]*(dl_th[0,*]+1))
        cl_th[3,*] = dl_th[3,*]*2*!DPI/(dl_th[0,*]*(dl_th[0,*]+1))

        spectrum = cl_th
    endif    





print, 'finished getting highel theory spectrum'

if keyword_set(stopit) then stop
return, spectrum
END

