pro setup_eighth_runlists

fields = ['ra3h30dec-60','ra21hdec-50','ra21hdec-60']


for i=0,2 do begin
    file='/home/cr/runlists/ltpairs/timeordered/setlist_coadd4X_timeordered_goodfiles_ppr_'+fields[i]+'.txt'
    
    readcol,file,col1,col2,col3,col4,format='a,a,a,a'
    
    crows = col1[0:*:8]
    
    ofile ='/home/cr/runlists/eighthlist_'+fields[i]+'.txt'
    
    n = n_elements(crows)
    openw,lun,/get_lun,ofile
    for j=0,n-1 do begin
        printf,lun,col1[4*j]
        printf,lun,col2[4*j]
        printf,lun,col3[4*j]
        printf,lun,col4[4*j]

    endfor
    free_lun,lun
endfor

end
