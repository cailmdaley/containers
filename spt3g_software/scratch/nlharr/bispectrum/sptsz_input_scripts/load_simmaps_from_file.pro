function load_simmaps_from_file, file_name, field
;must set if changed
nsims = 100
;;;;;;;;;;;;;;;;;;;;;;

fields_struct = get_field_struct()
fs = fields_struct[where(fields_struct.name eq field)]

print, "loading sim map with dims", fs.n1, fs.n2, nsims
out_array = read_binary(file_name, DATA_DIMS=[fs.n1, fs.n2, nsims], DATA_TYPE=4) ;data_type=4 float, 5 double

return, out_array

end
