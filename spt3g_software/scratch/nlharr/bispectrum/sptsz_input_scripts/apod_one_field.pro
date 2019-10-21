function apod_one_field, field, coadd_file, out_file, size_increase

print, 'field', field
print, 'coadd_file', coadd_file
print, 'out_file', out_file
print, 'size_increase', size_increase


;gets the info for the field
field_structs = get_field_struct()
fs = field_structs[where( field_structs.name eq field)]

;sets the parameters the function needs
reso_arcmin = 0.5
proj = 5

;help, fs, /st
;help, field_structs, /st
ra0 = fs.ra0
dec0 = fs.dec0

;reads the coadd_file
;fts = read_spt_fits(coadd_file)
;fts = spt_expand_smallmap(fts)
;wght = fts.weight.map
restore, coadd_file
wght = coadd.weight

create_apodization_masks, wght, '', reso_arcmin=reso_arcmin, mapproj = proj, radec0 = [ra0, dec0], output_struct=out_struct

mask_small = out_struct.apod_mask
ap_size = size(mask_small, /dimensions)

sq_side_len = long(ceil(max(size(mask_small, /dimensions)) * size_increase))
sq_side_len += sq_side_len mod 2 ;makes it so the mask is even and square

mask = make_array( sq_side_len, sq_side_len, value = 0.0)
mask[0:ap_size[0]-1,0:ap_size[1]-1] = mask_small

if n_elements(out_file) gt 0 then begin
    if out_file ne '' then begin
        save, mask, filename = out_file
    endif
endif

return, mask
end

