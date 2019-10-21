PRO SIM_maxsz,isim

sisim=strtrim(string(isim),2)
; Generate roughly R11-like simulated cluster catalogs
root = '/data/cr/paramfits/cosmomc.cluster/sims/'

restore, root+'params_highres.sav'
restore, root+'csm.sav'
restore, root+'sr.sav'
;restore, root+'cat.sav'
field_names = ['all']

field_areas = [2500.]

sr.asz*=1.2
field_depths = [1.]

if keyword_set(unidepth) then field_depths[*]=1d0
if keyword_set(hiuni) then field_depths[*]=1.2
field_asz = SR.asz * field_depths
nfld = n_elements(field_names)
; generate each of the catalogs
for i=0,nfld-1 do begin
  PARAMS.area = field_areas[i]
  SR.asz = field_asz[i]
  print,i,field_areas[i],field_asz[i]
  calcosmo_generate_mf, PARAMS, CSM, MF
  calcosmo_generate_catalog, CSM, PARAMS, SR, MF, CAT, /make_profile, loud=0
  out_struct = { field_name:field_names[i], $
                 field_area:field_areas[i], $
                 field_depth:field_depths[i], $
                 field_asz:field_asz[i], $
                 cat:CAT}
  dum = EXECUTE("os"+strtrim(i+1,2)+" = out_struct")
endfor

ns = intarr(1)
ns[0]=n_elements(os1.cat.x_vec)

n = total(ns)

n2 = n_elements(os1.CAT.txr_vec[0,*])

cat = {cat_name:'sim_'+sisim,ncluster:ns,$
       areas:field_areas,scalefac:field_depths,$
       falsefiles:['/home/cr/data/catalog/falserate/ra3h30dec-60_fr.txt'],$
       sz_vec:fltarr(n),$
       x_vec:fltarr(n),$
       logx_err_vec:fltarr(n),$
       z_vec:fltarr(n),$
       z_err_vec:fltarr(n),$
       da_vec:fltarr(n),$
       m_vec:fltarr(n),$
       mgr_vec:fltarr(n,n2),$
       txr_vec:fltarr(n,n2),$
       r_mg:fltarr(n,n2),$
       r_tx:fltarr(n,n2)}
icurr=0    


for i=0,nfld-1 do begin
   si = strtrim(string(i+1),2)
   dum = EXECUTE("tmp = os"+si+".cat")
   ncurr=ns[i]
   
   cat.sz_vec[icurr:icurr+ncurr-1]=tmp.sz_vec
   cat.x_vec[icurr:icurr+ncurr-1]=tmp.x_vec
   cat.logx_err_vec[icurr:icurr+ncurr-1]=tmp.logx_err_vec
   cat.z_vec[icurr:icurr+ncurr-1]=tmp.z_vec
   cat.z_err_vec[icurr:icurr+ncurr-1]=tmp.z_err_vec
   cat.m_vec[icurr:icurr+ncurr-1]=tmp.m_vec
   cat.da_vec[icurr:icurr+ncurr-1]=tmp.da_vec
   cat.mgr_vec[icurr:icurr+ncurr-1,*]=tmp.mgr_vec
   cat.txr_vec[icurr:icurr+ncurr-1,*]=tmp.txr_vec
   cat.r_mg[icurr:icurr+ncurr-1,*]=tmp.r_mg
   cat.r_tx[icurr:icurr+ncurr-1,*]=tmp.r_tx
   icurr+=ncurr
endfor


base='simcatmaxsz'

calcosmo_catalog_to_text,cat,root+base+strtrim(string(isim),2)+'.txt'
catalog=cat


save, catalog,os1, $
      filename=root+base+strtrim(string(isim),2)+'.sav'
;stop
RETURN
END
