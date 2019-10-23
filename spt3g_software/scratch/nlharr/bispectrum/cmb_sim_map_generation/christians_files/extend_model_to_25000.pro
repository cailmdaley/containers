pro extend_model_to_25000,ifile,ofile,norm=norm


readcol,ifile,l,dl,format='l,d'
nl = n_elements(l)
ml = l[nl-1]
mlp = ml-1000
dlm=dl[nl-1]
dlp=interpol(dl,l,mlp)

slope = (dlm-dlp)/(ml-mlp)

newl=lindgen(25001)
newdl=dblarr(25001)
newdl[l]=dl
newdl[nl:25000] = dlm+slope*lindgen(25001-nl)
ind = where(newdl lt 0,nii)
if nii gt 0 then newdl(ind)=0.0

if keyword_set(norm) then begin
    newdl/=newdl[3000]
endif

openw,1,ofile
for i=0l,25000 do printf,1,newl[i],newdl[i],format='(i,e)'
close,1
end
