

pro make_all_fields_initial
out_dir = '/data/nlharr/bispec_initial_run/'
tmp_tag = 'initial'
cluster_mask_size = 3
big_ptsrc_mask_size = 3
small_ptsrc_mask_size = 2
map_size_increase = 1.0
cluster_mass_cut = 2.5e-14

;800 sq deg fields; 


field_lst = ['ra5h30dec-45', 'ra5h30dec-45', 'ra5h30dec-45', 'ra2h30dec-50', 'ra2h30dec-50', 'ra2h30dec-50', 'ra1hdec-60', 'ra1hdec-60', 'ra1hdec-60', 'ra4h10dec-50', 'ra4h10dec-50', 'ra4h10dec-50', 'ra21hdec-42.5', 'ra21hdec-42.5', 'ra21hdec-42.5', 'ra6h30dec-55', 'ra6h30dec-55', 'ra6h30dec-55', 'ra3h30dec-42.5', 'ra3h30dec-42.5', 'ra3h30dec-42.5', 'ra23hdec-62.5', 'ra23hdec-62.5', 'ra23hdec-62.5', 'ra23hdec-45', 'ra23hdec-45', 'ra23hdec-45', 'ra22h30dec-55', 'ra22h30dec-55', 'ra22h30dec-55', 'ra23h30dec-55', 'ra23h30dec-55', 'ra23h30dec-55', 'ra1hdec-42.5', 'ra1hdec-42.5', 'ra1hdec-42.5', 'ra21hdec-50', 'ra21hdec-50', 'ra21hdec-50', 'ra3h30dec-60', 'ra3h30dec-60', 'ra3h30dec-60', 'ra5h30dec-55', 'ra5h30dec-55', 'ra5h30dec-55', 'ra6hdec-62.5', 'ra6hdec-62.5', 'ra6hdec-62.5', 'ra6h30dec-45', 'ra6h30dec-45', 'ra6h30dec-45', 'ra21hdec-60', 'ra21hdec-60', 'ra21hdec-60', 'ra0h50dec-50', 'ra0h50dec-50', 'ra0h50dec-50']

band_lst = [220, 90, 150, 220, 150, 90, 150, 220, 90, 220, 150, 90, 220, 90, 150, 220, 90, 150, 150, 220, 90, 220, 150, 90, 90, 220, 150, 150, 90, 220, 150, 90, 220, 220, 90, 150, 150, 220, 90, 150, 90, 220, 90, 150, 220, 150, 220, 90, 220, 150, 90, 90, 150, 220, 220, 150, 90]


;field_lst = ['ra5h30dec-55']
;band_lst = [150]


len_fbs = size(field_lst, /dim)
len_fbs = len_fbs[0]


for ii=0, len_fbs-1 do begin
  make_complete_one_field, field_lst[ii], out_dir, tmp_tag, band_lst[ii], cluster_mask_size, big_ptsrc_mask_size, small_ptsrc_mask_size, map_size_increase, cluster_mass_cut, 1
endfor
end
