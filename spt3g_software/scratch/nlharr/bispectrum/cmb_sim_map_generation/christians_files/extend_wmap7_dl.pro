pro extend_wmap7_dl


readcol,'~cr/cmb_models/wmap7_lcdm_lensedCls.dat',l,dl

inds = where(l gt 7500)
lx = l(inds)
lndl = alog(dl(inds))
res = linfit(lx,lndl)

lp = max(l)+1+findgen(10000)
lndlp = res[0]+res[1]*lp
dlp = exp(lndlp)
ll = [l,lp]
dll = [dl,dlp]

openw,1,'~cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat'
nl = n_elements(ll)
for i=0l,nl-1 do printf,1,format='(i,e)',ll[i],dll[i]
close,1
;stop

end
