pro quick_guess_spt_full_bandpowers,es=es

if n_elements(es) ne 1 then es = ''
make_ml_model_curves_fixsz,'~/data/chains_r11/chains_20110529/chain_r11l_sptlowl9p2_hfi217_shawtsz_csfksz_hall_1.txt',3.5,2.0,0,ml_l,ml_dls,ml_cal

dls = reform(ml_dls[*,0,*]+ml_dls[*,1,*])
;cut to winf ells 50+
dls = dls[48:*,*]

wstub='~/r11/window_'
winfs = dblarr(10999,90)
lwin=50+findgen(10999)
for i=0,89 do begin
    a = read_ascii(wstub+strtrim(string(i+1),2))
    winfs[*,i]=a.field1[1,*]
endfor
;cut to dls ells: <=10000
winfs=winfs[0:9950,*]

;now produce bandpowers
bdls = dblarr(90)
for i=0,89 do begin
    j = i / 15
    bdls[i] = total(dls[*,j]*winfs[*,i])
endfor

;now for errors...
a = read_ascii('~cr/r11_newcalbm/covariance_90_90x150_90x220_150_150x220_220.txt')
cov = reform(a.field1,90,90)
cov /= 3.0


errvec = realize_noise_cov(cov)

;double check this before running
;stop
bdls0=bdls
bdls += errvec
print,moment(errvec)
;stop

;now write to disk
odir = '~cr/futurespt'+es+'/'
spawn,'mkdir '+odir
sfile='spectrum_90_90x150_90x220_150_150x220_220.txt'
cfile='covariance_90_90x150_90x220_150_150x220_220.txt'
openw,lun,/get_lun,odir+sfile
for i=0,89 do $
  printf,lun,format='(i,e)',i mod 15, bdls[i]
close,lun
openw,lun,odir+cfile
for i=0,89 do for j=0,89 do $
  printf,lun,format='(e)',cov(j,i)
free_lun,lun
stop
end
