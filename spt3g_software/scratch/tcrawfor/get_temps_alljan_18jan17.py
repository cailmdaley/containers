import glob
import get_temperatures

arcdir = '/spt_data/arc/'
arcfiles = np.asarray(glob.glob(arcdir+'201801*.dat'))
arcfiles.sort()
arcfiles = arcfiles[350:] # only use from about 1/5
outdir = '/home/tcrawfor/temps_05jan17_to_17jan17/'

for file1 in arcfiles:
    print(file1)
    print("")
    outfile = file1.replace(arcdir,outdir)
    outfile = outfile.replace('dat','pkl')
    try:
        dtemp = get_temperatures.get_T(file1,outfile)
    except:
        print("Couldn't get temps from file "+file1)
    print("")
