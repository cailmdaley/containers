pro write_suman_tsz

ll=dindgen(25000)+2
a = read_ascii('~cr/cmb_models/suman_tsz_z4_20131220.txt')
tsz=ll*0.0
tszdl=[0,reform(a.field1[1,*]),0]
tszl=[0,reform(a.field1[0,*]),25001]

tsz = interpol(tszdl,tszl,ll)
tsz[0:9]=0
tsz/=tsz[3000]

openw,1,'~cr/cmb_models/dl_tsz_suman_20131220_norm1_l25000.txt'
for i=0l,24999 do begin
printf,1,ll[i],tsz[i],format='(i,e)'
endfor

close,1
plot,ll,tsz
;stop
end
