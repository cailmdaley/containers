

function get_field_obs_lst
flst_90 = ['ra0h50dec-50', 'ra1hdec-42.5', 'ra1hdec-60', 'ra21hdec-42.5', 'ra21hdec-50_elscan', 'ra21hdec-50', 'ra21hdec-60', 'ra22h30dec-55', 'ra23h30dec-55_2010', 'ra23hdec-45', 'ra23hdec-62.5', 'ra2h30dec-50', 'ra3h30dec-42.5', 'ra3h30dec-60', 'ra4h10dec-50', 'ra5h30dec-45', 'ra5h30dec-55_2011', 'ra6h30dec-45', 'ra6h30dec-55', 'ra6hdec-62.5']
flst_150 = ['ra1hdec-60', 'ra21hdec-50_elscan', 'ra23h30dec-55_2008', 'ra23hdec-45', 'ra3h30dec-42.5', 'ra5h30dec-45', 'ra6h30dec-45', 'ra0h50dec-50', 'ra21hdec-42.5', 'ra21hdec-60', 'ra23h30dec-55_2010', 'ra23hdec-62.5', 'ra3h30dec-60', 'ra5h30dec-55_2008', 'ra6h30dec-55', 'ra1hdec-42.5', 'ra21hdec-50', 'ra22h30dec-55', 'ra23h30dec-55_elscan', 'ra2h30dec-50', 'ra4h10dec-50', 'ra5h30dec-55_2011', 'ra6hdec-62.5']
flst_220 = ['ra0h50dec-50', 'ra1hdec-42.5', 'ra1hdec-60', 'ra21hdec-42.5', 'ra21hdec-50_elscan', 'ra21hdec-50', 'ra21hdec-60', 'ra22h30dec-55', 'ra23h30dec-55_2008', 'ra23h30dec-55_2010', 'ra23h30dec-55_elscan', 'ra23hdec-45', 'ra23hdec-62.5', 'ra2h30dec-50', 'ra3h30dec-42.5', 'ra3h30dec-60', 'ra4h10dec-50', 'ra5h30dec-45', 'ra5h30dec-55_2008', 'ra5h30dec-55_2011', 'ra6h30dec-45', 'ra6h30dec-55', 'ra6hdec-62.5']

return, {f90:flst_90, f150:flst_150, f220:flst_220}

end
