pro write_sim09_ksz

ll=dindgen(25000)+2
a = read_ascii('~cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt')
ksz=ll*0.0
ksz[0:16300]=a.field1[1,2:16302]

;smoothly go to zero                                                            
ksz[16301:24999] = ksz[16300] - (ksz[16300]/8699.)*(1+findgen(8699))
;l=2,25000

openw,1,'~cr/cmb_models/dl_ksz_sehgal_sim09.txt'
for i=0l,24999 do begin
printf,1,ll[i],ksz[i],format='(i,e)'

endfor

close,1


end
