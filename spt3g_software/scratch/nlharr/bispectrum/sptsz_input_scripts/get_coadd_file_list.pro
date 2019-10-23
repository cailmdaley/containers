

;example field name is: ra0h50dec-50

function get_coadd_file_list, field_name, band
band_str = strtrim(string(band), 2)
return, file_search('/data23/lizinvt/highel2500/coadds/combmapsx2_jittertimeordered_20130616_good/'+field_name+'/map_'+field_name+'_'+band_str+'*fits')

rlst = get_run_list(field_name)

;return, rlst
;end
;; reads the runlist into an array
OPENR, lun, rlst, /GET_LUN
; Read one line at a time, saving the result into array

line = ''
WHILE (~ EOF(lun)) DO BEGIN
  READF, lun, line
  if n_elements(rlarray) eq 0 then rlarray=[line] else rlarray=[rlarray,line]
ENDWHILE
; Close the file and free the file unit
FREE_LUN, lun


out_files = '/data23/lizinvt/highel2500/data_single_maps/'+field_name+'/maps/map_'+field_name+'_'+strtrim(string(band),2)+'_'+rlarray+'.fits'
;/data23/lizinvt/highel2500/data_single_maps/ra0h50dec-50/maps/map_ra0h50dec-50_150_20100618_010104.fits
return, out_files



;band_str = strtrim(string(band), 2)
;return, file_search('/data23/lizinvt/highel2500/coadds/combmapsx2_jittertimeordered_20130616_good/'+field_name+'/map_'+field_name+'_'+band_str+'*fits')
end
