function spt_coadd_dir,dir,tag,v1=v1,v2=v2,nowt=nowt
; a very brief code to coadd all sav files in a directory
; includes almost no error checking and will crash if there are no files.

list = file_search(dir+tag+'*')
;stop
nlist = n_elements(list)
;list=list(nlist/2:*)
;nlist=nlist-nlist/2

if nlist le 1 then return,{ncoadd:nlist}

map = 0
nobs=0
if keyword_set(nowt) then begin
    
    if keyword_set(v1)+keyword_set(v2) eq 0 then begin
        m = read_spt_fits(list[0])
        if find_matching_tag(m,'PIXELS') eq 'PIXELS' then begin
            map = dblarr(m.mapinfo.nsidex,m.mapinfo.nsidey)
            nobs = map
            for i=0,nlist-1 do begin
;                print,i,nlist
                m = read_spt_fits(list[i])
                pixels = convert_pixel_array(m)
                map(pixels) += m.map.map
                inds =where(m.weight.map gt 0)
                nobs(pixels(inds)) += 1.
            endfor

        endif else $
          for i=0,nlist-1 do begin
            m = read_spt_fits(list[i])
            map = map+m.map.map
            if i eq 0 then nobs = map*0.0
            inds =where(m.weight.map gt 0)
            nobs(inds) += 1.
        endfor
    endif else $
      if ~ keyword_set(v1) then begin
        for i=0,nlist-1 do begin
            restore,list[i]
            map = map+mapstruct.map
            if i eq 0 then nobs = map*0.0
            inds =where(mapstruct.weight gt 0)
            nobs(inds) += 1.
        endfor
    endif else begin
        for i=0,nlist-1 do begin
            restore,list[i]
            map = map+smap.map
            if i eq 0 then nobs = map*0.0
            inds =where(smap.weight gt 0)
            nobs(inds) += 1.
        endfor
    endelse
endif else begin
    if keyword_set(v1)+keyword_set(v2) eq 0 then begin
        m = read_spt_fits(list[0])
        if find_matching_tag(m,'PIXELS') eq 'PIXELS' then begin
            map = dblarr(m.mapinfo.nsidex,m.mapinfo.nsidey)
            nobs = map
            for i=0,nlist-1 do begin
;                print,i,nlist
                m = read_spt_fits(list[i])
                pixels = convert_pixel_array(m)
                map(pixels) += m.map.map*m.weight.map
                nobs(pixels) += m.weight.map
            endfor
        endif else $
          for i=0,nlist-1 do begin
            m = read_spt_fits(list[i])
            map = map+m.map.map*m.weight.map
            nobs += m.weight.map
        endfor
    endif else $
      if ~ keyword_set(v1) then begin
        for i=0,nlist-1 do begin
            restore,list[i]
            map = map+mapstruct.map*mapstruct.weight
            nobs = nobs+mapstruct.weight
        endfor
    endif else begin
        for i=0,nlist-1 do begin
            restore,list[i]
            map = map+smap.map*smap.weight
            nobs = nobs+smap.weight
        endfor
    endelse
    
endelse



wh = where(nobs gt 0)
map[wh]=map[wh]/nobs[wh]

return,{map:map,weight:nobs,ncoadd:nlist,runlist:list}
end
