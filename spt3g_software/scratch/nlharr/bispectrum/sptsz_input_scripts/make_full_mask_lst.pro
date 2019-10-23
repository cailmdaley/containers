
;use a cutoff of 2.5e14


;example call:
; make_full_mask_lst, "/home/nlharr/Bispec_2500/pt_src_lists/ptsrc_config_ra0h50dec-50.txt", "ra0h50dec-50", "/home/nlharr/Bispec_2500/pt_src_lists/spt_clust_allz_list_07mar2014.txt", 2.5e-14, full_ra, full_dec, full_desc

pro make_full_mask_lst, ptsrc_file, field, cluster_file, cluster_mass_cut, full_ra, full_dec, full_desc
;I'm gonna decide right now that there is a descriptor of the point src type

;0 = cluster
;1 = big ptsrc
;2 = small ptsrc


reso_arcmin = 0.5
proj = 5


;loads the field information
field_struct = get_field_struct()
fs = field_struct[where( field_struct.name eq field)]
ra0 = fs.ra0
dec0 = fs.dec0
n2 = fs.n2
n1 = fs.n1

radec0 = [ra0, dec0]
npixels = [n1,n2]

;read the point source file
print, 'reading ptsrc file ', ptsrc_file
readcol, ptsrc_file, index, ptsrc_ra, ptsrc_dec, ptsrc_radius, SKIPLINE=3, COMMENT='#', format='IDDD'

;0 = cluster
;1 = big ptsrc
;2 = small ptsrc

ptsrc_descriptor = long(ptsrc_radius) * 0 + 2
ptsrc_descriptor[where(ptsrc_radius gt 0.05)] = 1

print, ptsrc_radius


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;load the clusters
;db = check_database(/official)
;;;;filter them based off of validity criteria
;val_clusters = db[where(db.spt_confidence gt 5)]
;val_cluster_ras = val_clusters.spt_ra
;val_cluster_decs = val_clusters.spt_dec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

parse_cluster_list_file, cluster_file, cluster_mass_cut, val_cluster_ras, val_cluster_decs

;filters them based off of which ones are in the field
ang2pix_proj5, val_cluster_ras, val_cluster_decs, npixels, radec0, reso_arcmin, out_pix

valval_cluster_inds = where(out_pix gt -1)


if size(valval_cluster_inds, /N_DIMENSIONS) ne 0 then begin
    fin_ras = val_cluster_ras[valval_cluster_inds]
    fin_decs = val_cluster_decs[valval_cluster_inds]
    cluster_descriptor = long(fin_ras) * 0
    full_ra  = [ptsrc_ra, fin_ras]
    full_dec = [ptsrc_dec, fin_decs]
    full_desc = [ptsrc_descriptor, cluster_descriptor]
endif else begin
    full_ra = ptsrc_ra
    full_dec = ptsrc_dec
    full_desc = ptsrc_descriptor
endelse





end
