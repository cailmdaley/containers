

pro parse_cluster_list_file, file_name, cluster_mass_cut, ras, decs
readcol, file_name, cluster_name, ra, dec, z, sn, m500, format="ADDDDD"

cluster_inds = where(m500 gt cluster_mass_cut)
ras = ra[cluster_inds]
decs = dec[cluster_inds]

end
