import glob, filecmp, copy

def uniquifyList(seq):
    '''
    Returns a list with all the elements of seq where all the duplicate elements are removed
    '''
    # order preserving
    checked = []
    for e in seq:
        if e not in checked:
            checked.append(e)
    return checked


fields = ['ra0h50dec-50', 'ra1hdec-42.5', 'ra1hdec-60', 'ra21hdec-42.5', 'ra21hdec-50', 'ra21hdec-60', 'ra22h30dec-55', 'ra23h30dec-55', 'ra23hdec-45', 'ra23hdec-62.5', 'ra2h30dec-50', 'ra3h30dec-42.5', 'ra3h30dec-60', 'ra4h10dec-50', 'ra5h30dec-45', 'ra5h30dec-55', 'ra6h30dec-45', 'ra6h30dec-55', 'ra6hdec-62.5']
freqs = ['90','150','220']

check_fmt = '/data30/nlharr/Bispectrum2500D/initial_run_output/out_file_${field}_${freq}GHz_*.hdf5'

unique_keys = []
for field in fields:
    for freq in freqs:
        unique_keys.append([field,freq])
    
num_files = []        
for uk in unique_keys:
    glob_str = check_fmt.replace('${field}', uk[0]).replace('${freq}', uk[1])
    num_files.append(len(glob.glob(glob_str)))


for i in xrange(len(unique_keys)):
    print unique_keys[i], num_files[i]

        
full_glob_str = check_fmt.replace('${field}', '*').replace('${freq}', '*')
flst = glob.glob(full_glob_str)


print 'checking that none are the same'
for i, f in enumerate(flst):
    for j, g in enumerate(flst[i+1:]):
        assert(not filecmp.cmp(f,g))




