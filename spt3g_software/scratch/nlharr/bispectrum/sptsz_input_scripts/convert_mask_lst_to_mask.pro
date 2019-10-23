


;the mask size is in units of pixels away from the center we mask,
;where the size of the pixels we mask is 1 + 2*size.  So size = 2
;means masking a five pixel square

;code to test this:

;make_full_mask_lst, "/home/nlharr/Bispec_2500/pt_src_lists/ptsrc_config_ra0h50dec-50.txt", "ra0h50dec-50", "/home/nlharr/Bispec_2500/pt_src_lists/spt_clust_allz_list_07mar2014.txt", 2.5e-14, full_ra, full_dec, full_desc
;mask = convert_mask_lst_to_mask(full_ra, full_dec, full_desc, "ra0h50dec-50", 5, 5, 3)


function convert_mask_lst_to_mask, ras, decs, descs, field, cluster_mask_size, big_ptsrc_mask_size, small_ptsrc_mask_size

reso_arcmin = 0.5


;uses the size of the apod mask to generate the ptsrc mask


;gets the ra and dec vals of the field
field_structs = get_field_struct()
fs = field_structs[where( field_structs.name eq field)]

ap_size = size(apod_mask)
d1 = fs.n1
d2 = fs.n2
mask = intarr( d1,d2) * 0 ;just making sure it's zero because I am bad at idl


;
ra0 = fs.ra0
dec0 = fs.dec0

;iterates through input list and masks based off of input parameters

n_ptsrc_tmp = size(ras)

n_ptsrc = n_ptsrc_tmp[1]

;first we get the indices
radec0 = [ra0, dec0]
npixels = [fs.n1, fs.n2]
ang2pix_proj5, ras, decs,npixels, radec0, reso_arcmin, xpix=xpix, ypix=ypix


FOR ii = 0, n_ptsrc-1 DO BEGIN
    if descs[ii] eq 0 then pt_mask_size = cluster_mask_size
    if descs[ii] eq 1 then pt_mask_size = big_ptsrc_mask_size
    if descs[ii] eq 2 then pt_mask_size = small_ptsrc_mask_size
    ;masks based off of size
    for jj = -1*pt_mask_size, pt_mask_size do begin
        for kk = -1*pt_mask_size, pt_mask_size do begin
            mask[jj+xpix,kk+ypix] = 1
        endfor
    endfor
ENDFOR


return, mask

end
