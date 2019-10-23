pro breakup_run4_2009_sims
nper=100

year='2008'
freqs=['150','220']
fields=['ra23h30dec-55','ra5h30dec-55']
n1=3228
n2=2988
n3=600
for i=0,n_elements(freqs)-1 do for j=0,n_elements(fields)-1 do begin
    field=fields[j]
    freq=freqs[i]
    filein = '/data11/cr/r11sims/hrsim_spt_run4_2009_map'+year+'_'+freq+'_gauss.dat'
    fileout = '/data11/cr/r11sims/rotfhrsim_spt_run4_2009_map'+year+'_'+freq+'_gauss_'+field+'.dat'
    
;    rotate_partial_sim_map,filein,fileout,n1,n2,n3,on3=nper*j,un3=100
endfor

year='2009'
freqs=['90','150','220']
fields=['ra21hdec-50','ra21hdec-60','ra3h30dec-60']
n1=6612
n2=3512
n3=900
for i=0,n_elements(freqs)-1 do for j=0,n_elements(fields)-1 do begin
    field=fields[j]
    freq=freqs[i]
    filein = '/data11/cr/r11sims/hrsim_spt_run4_2009_map'+year+'_'+freq+'_gauss.dat'
    fileout = '/data11/cr/r11sims/rotfhrsim_spt_run4_2009_map'+year+'_'+freq+'_gauss_'+field+'.dat'
    
    rotate_partial_sim_map,filein,fileout,n1,n2,n3,on3=nper*j,un3=100
endfor

end
