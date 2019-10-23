# how to even out sub-field noise for 2019 (starting in August)

# which band to do the evening in
band2use = '90'
#band2use = '150'

# noise numbers as of 7/29 from Wei
bands = ['90','150','220']
fields = ['ra0hdec-44.75','ra0hdec-52.25','ra0hdec-59.75','ra0hdec-67.25']
noise_dict = {}
for band in bands:
    noise_dict[band] = {}
    for field in fields:
        noise_dict[band][field] = {}
    noise_dict[band]['ra0hdec-44.75']['nmaps'] = 210
    noise_dict[band]['ra0hdec-52.25']['nmaps'] = 196
    noise_dict[band]['ra0hdec-59.75']['nmaps'] = 195
    noise_dict[band]['ra0hdec-67.25']['nmaps'] = 169
nmaps_done = 0.0
for field in fields:
    nmaps_done += np.float(noise_dict['90'][field]['nmaps'])
noise_dict['90']['ra0hdec-44.75']['noise'] = 11.9
noise_dict['90']['ra0hdec-52.25']['noise'] = 11.5
noise_dict['90']['ra0hdec-59.75']['noise'] = 11.2
noise_dict['90']['ra0hdec-67.25']['noise'] = 10.7
noise_dict['150']['ra0hdec-44.75']['noise'] = 9.9
noise_dict['150']['ra0hdec-52.25']['noise'] = 9.1
noise_dict['150']['ra0hdec-59.75']['noise'] = 8.9
noise_dict['150']['ra0hdec-67.25']['noise'] = 8.3
noise_dict['220']['ra0hdec-44.75']['noise'] = 36.6
noise_dict['220']['ra0hdec-52.25']['noise'] = 34.5
noise_dict['220']['ra0hdec-59.75']['noise'] = 32.0
noise_dict['220']['ra0hdec-67.25']['noise'] = 30.6
for band in bands:
    for field in fields:
        noise_dict[band][field]['var_single_obs'] = noise_dict[band][field]['noise']**2*noise_dict[band][field]['nmaps']

# figure out how many observations we are getting per week
#  ok, an uninterrupted "week" of ABBABBCD seems to have taken this long:
week_start = '190707 20:51:59'
week_end = '190714 11:12:57'
week_duration_days = core.G3Time(week_end).mjd - core.G3Time(week_start).mjd
# and we get 45 field obs in that time
nmaps_per_day = 45./week_duration_days


# so how many observations total are left in the season?
##  first assume observing through 10/31
#ndays_left = core.G3Time('20191031_000000').mjd - core.G3Time('20190731_000000').mjd
#nmaps_left = ndays_left*nmaps_per_day
#nmaps_tot = nmaps_done + nmaps_left
##  now assume observing through 11/31
#ndays_left = core.G3Time('20191131_000000').mjd - core.G3Time('20190731_000000').mjd
#nmaps_left = ndays_left*nmaps_per_day
#nmaps_tot = nmaps_done + nmaps_left
#  Brad's best guess: November 13
ndays_left = core.G3Time('20191113_000000').mjd - core.G3Time('20190731_000000').mjd
nmaps_left = ndays_left*nmaps_per_day
nmaps_tot = nmaps_done + nmaps_left

# we want the total number of maps in each field to be proportional to the single-obs variance in that field
var_tot = 0.
nmaps_tot_desired = {}
nmaps_to_do = {}
for field in fields:
    var_tot += noise_dict[band2use][field]['var_single_obs']
for field in fields:
    nmaps_tot_desired[field] = noise_dict[band2use][field]['var_single_obs']*nmaps_tot/var_tot
    nmaps_to_do[field] = np.round(nmaps_tot_desired[field]-noise_dict[band2use][field]['nmaps'])



