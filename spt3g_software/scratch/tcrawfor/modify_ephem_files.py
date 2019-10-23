# script to take ephem output from HORIZONS and get into format GCP
# seems to want

dir = '/home/tcrawfor/temp_ephem_2017/'
#target_list = ['sun','moon','mercury','venus','mars','jupiter','saturn','uranus','neptune','pluto']
target_list = ['mercury']

for target in target_list:
    f1 = open(dir + target + '_unmod.ephem','r')
    f2 = open(dir + target + '.ephem','w')
    f2.write('#                     SPT Ephemeris for ' + target  + '\n')
    f2.write('#\n')
    f2.write('# Note that TT ~= TAI + 32.184.\n')
    f2.write('# See slaDtt.c for how to convert TT to UTC.\n')
    f2.write('#\n')
    f2.write('# Also note that MJD = JD - 2400000.5\n')
    f2.write('#\n')
    f2.write('# MJD (TT)   Right Ascen    Declination     Distance (au)       Calendar (TT)\n')
    f2.write('#---------  -------------  -------------   ----------------    ------------------------\n')
    for line in f1.readlines():
        fline = filter(None,line.split(' '))
        mjdstr = '%10.4f' % (np.float(fline[2]) - 2400000.5)
#        rastr = fline[3] + ':' + fline[4] + ':' + '%02.4f' % (np.float(fline[5]))
        rahstr = fline[3]
        ramstr = '%02d' % np.int(fline[4])
        rasstr = '%07.4f' % (np.float(fline[5]))
        rastr = rahstr + ':' + ramstr + ':' + rasstr

        dechstr = fline[6]
        decmstr = '%02d' % np.int(fline[7])
        decsstr = '%06.3f' % (np.float(fline[8]))
        decstr = dechstr + ':' + decmstr + ':' + decsstr

        austr = fline[9]
        calstr = fline[0] + ' ' + fline[1]
        f2.write(mjdstr + '  ' + rastr + '  ' + decstr + '   ' + austr + ' #  ' + calstr + '\n')

    f2.close()
    f1.close()
