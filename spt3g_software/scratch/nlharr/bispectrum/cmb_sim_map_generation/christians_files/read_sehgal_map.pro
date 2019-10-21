function read_sehgal_map,file
;this pro fails to work!
@boom_com
getindata, [file]
;stop
;units Jy/steradion
; I want mJy
;nside=8192
nall=long64(12)*long64(nside)*nside
pixarea=4.*!pi/nall
signal *=pixarea *1000.; now in mJy
;inds = where(signal gt 6.7,ni)
return,signal

end
