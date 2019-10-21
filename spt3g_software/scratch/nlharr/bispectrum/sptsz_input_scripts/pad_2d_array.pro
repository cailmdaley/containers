



function pad_2d_array, array, d0, d1
  out_arr = make_array( d0, d1, value = 0.0)
  arr_size = size(array, /dimension)
  ;print, 'arr_size', arr_size
  out_arr[0:arr_size[0]-1,0:arr_size[1]-1] = array
  return, out_arr
end
