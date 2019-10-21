import numpy as np
from spt3g.util import diff_scan_parser
from spt3g.util import diff_scan_maker

gcpdir = '/home/tcrawfor/code/GCP/'
schdir = gcpdir + 'config/sch/'
sched3 = schdir + 'scan_1500d_field_elevation_3_13Feb2018.sch'

el_max = 70.
el_min = 42.
ra_max = 50.
ra_min = -50.
el_buffer = 1.
xel_buffer = 1.
el_throw = el_max - el_min + 2.*el_buffer

nsched = 4
el_throw_sched = el_throw/nsched
el_center = np.arange(nsched)*el_throw_sched + el_throw_sched/2. + el_min - el_buffer
sec = np.asarray([str(np.int(np.round(elcen))) for elcen in el_center])
d2r = np.pi/180.
cec = np.cos(el_center*d2r)

# let's try the assumption that every schedule will have 60 el steps and we will do each twice in a fridge cycle
nelstep = 60
el_step_arcmin = el_throw_sched*60./nelstep
sesa = ("%3.1f" % el_step_arcmin).replace('.','p')
ndither = np.int(np.round(el_step_arcmin/0.5))
vscan = 1.9/cec*cec[3]
svscan = np.asarray([("%4.2f" % vs).replace('.','p') for vs in vscan])
amax = 0.38
samax = '0p38'
az_throw = ra_max - ra_min + np.asarray([2.*xel_buffer/tcec for tcec in cec])
saz_throw = np.asarray(["%3d" % np.int(np.round(ath)) for ath in az_throw])
scan_names = np.asarray(['daz-'+sat+'d-'+svs+'ds-'+samax+'dss-'+sesa+'am' for svs,sat in zip(svscan,saz_throw)])
offsets = []
scan_times = []
for ath, vs, sn in zip(az_throw,vscan,scan_names):
#    offsets.append(diff_scan_maker.scanConstJerk(ath, vs, amax, el_step_arcmin/60., filename=sn))
    scan = diff_scan_parser.Scan.fromFilename(sn+'.scan')
#    scan_times.append(np.sum(scan.flags)/100.)
    scan_times.append(float(len(scan.az))/100.)
scan_times = np.asarray(scan_times)

print ""
print "If you do each of these schedules twice, it will take %4.1f hours." % np.sum(scan_times*nelstep*2./3600.)
print ""

for i in np.arange(3):
    print "To write the schedule scan_1500d_field_elevation_0_13Feb2018.sch (assuming you are using scan_1500d_field_elevation_3_13Feb2018.sch as a template), you will need to replace the following lines:"
    print "  1) Replace all instances of 'scan_1500d_field_elevation_3_13Feb2018.sch' with 'scan_1500d_field_elevation_"+str(i)+"_13Feb2018.sch' (or just 'elevation_3' with 'elevation_"+str(i)+"')."
    print "  2) Replace all instances of '67.25d' with '"+("%5.2f" % el_center[i])+"d.'"
    print "  3) Replace all instances of 'el67' with 'el"+sec[i]+".'"
    if i < 2:
        print "  4) Comment out all lines referring to mat5a, and comment back in all lines referring to rcw38."
    print "  5) Replace the string '{daz-105d-1p9ds-0p38dss-7p5am, 60, -65.265, -3.75, 0.00833333, 15}' with the string '{"+scan_names[i]+", 60, "+("%7.3f" % offsets[i])+", -3.75, 0.00833333, 15}.'"
    print ""



#replace all #    scan_1500d_field_elevation_3_13Feb2018.sch
# This one does the highest elevation (centered at el=67.25d).
# centered at ra=0, dec=-67.25d.
# command scan_1500d_field_el67_2018(FieldObservation obs, Integer iEl_dither) {
#  log "In command scan_1500d_field_el67_2018 iEl_dither=", $iEl_dither 
# set_up_radio_cold_2018_el67p0
# set_up_radio_warm_2018_el67p0
# switch out mat5a/rcw38 based on elevation
#  scan_1500d_field_el67_2018 , $fastScan, $iEl_dither
#FieldObservation fastScan = {daz-105d-1p9ds-0p38dss-7p5am, 60, -65.265, -3.75, 0.00833333, 15}
