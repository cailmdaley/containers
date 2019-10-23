'''
This file contains functions related to the telescope pointing 
while the EHT (Event Horizon Telescope) is on SPT.  The most useful
function is get_thermolin_correction_now, which calculates small 
corrections to the collimation terms to be used in the online 
pointing model, and outputs the relevant command to be pasted 
into pointing.init.

Unless otherwise noted, all units are degrees.

created 08 Jan 2015, RK
'''

import numpy as np
np.random.seed(0)
import pdb
import pickle
import os
from spt3g import core, std_processing, gcp
#from sptpol_software.util import time
prefit_model_path = '/home/rkeisler/eht_pointing/'


def test(n_days_ago=100, base_col_x=0.3, base_col_y=-0.1):
    '''
    This function is just used for testing.
    Get the thermolin data from N_DAYS_AGO, 
    calculate the corrections to the collimation terms, and,
    assuming some (made-up) baseline collimation terms, print 
    out the new collimation terms (=baseline ones + correction).
    '''

    start_time = core.G3Time().Now() - n_days_ago*core.G3Units.days
    get_thermolin_correction(0.3, -0.1, start_time.Summary())

def get_thermolin_correction_now(base_col_x, base_col_y, 
                                 n_minutes_ago=5):
    '''
    This function calculates the correction to the collimation
    terms in the online model based on thermometry and linear
    sensor data ("thermolin" data) from "now", i.e. some number of 
    minutes ago.  The new "collimate_fixed" command is printed to screen.
    This function is basically just a wrapper to 
    GET_THERMOLIN_CORRECTION.

    Ideally, a new archive file was just written, so you can get 
    very recent (5 minute old) thermolin data.  The quality of the 
    correction will get worse as the thermolin data gets older, and I 
    wouldn't even try to use this correction if the thermolin data is
    older than 30 minutes.


    INPUTS

        BASE_COL_X - the baseline collimation in the X direction (i.e.
        the cross-elevation collimation).  The new collimation is 
        this one plus the thermolin-based correction. [degrees]

        BASE_COL_Y - the baseline collimation in the Y direction (i.e.
        the elevation collimation).  The new collimation is this one
        plus the thermolin-based correction. [degrees]


    OPTIONAL INPUTS

        N_MINUTES_AGO - We try to grab thermolin data from this many 
        minutes ago.  The fresher the better, and I wouldn't even 
        try to use this correction if the thermolin data is >30 min old.

    made 10 Jan 2015, RK
    '''
    
    start_time = core.G3Time().Now() - n_minutes_ago*core.G3Units.minutes
    get_thermolin_correction(base_col_x, base_col_y, start_time.Summary())

def get_thermolin_correction(base_col_x, base_col_y, start_time_):
    '''
    This function calculates the correction to the collimation
    terms in the online model based on thermometry and linear
    sensor data ("thermolin" data) from "now", i.e. some number of 
    minutes ago.  The new "collimate_fixed" command is printed to screen.

    Ideally, a new archive file was just written, so you can get 
    very recent (5 minute old) thermolin data.  The quality of the 
    correction will get worse as the thermolin data gets older, and I 
    wouldn't even try to use this correction if the thermolin data is
    older than 30 minutes.

    INPUTS

        BASE_COL_X - the baseline collimation in the X direction (i.e.
        the cross-elevation collimation).  The new collimation is 
        this one plus the thermolin-based correction. [degrees]

        BASE_COL_Y - the baseline collimation in the Y direction (i.e.
        the elevation collimation).  The new collimation is this one
        plus the thermolin-based correction. [degrees]

        START_TIME - a string (like '01-Jan-2015:00:00:00') or a 
            datetime object denoting the time you want thermolin 
            data from.  We then try to grab thermolin data from 
            a 10-second window starting on that start time.

    made 10 Jan 2015, RK.
    '''

    #from datetime import datetime, timedelta
    ## If START_TIME is a string, convert to a datetime object.
    #if isinstance(start_time_, str):
    #    fmt = '%d-%b-%Y:%H:%M:%S'
    #    start_time = datetime.strptime(start_time_, fmt)
    #else: start_time = start_time_
    start_time = core.G3Time(start_time_)
    stop_time = start_time + 10.*core.G3Units.seconds

    thermolin = read_thermolin_data_from_archive_3G(start_time.Summary(), stop_time.Summary())

    # Load the pre-fit models .
    rf_x = pickle.load(open(os.path.join(prefit_model_path,'thermolin_model_xoff.pkl'),'r'))
    rf_y = pickle.load(open(os.path.join(prefit_model_path,'thermolin_model_yoff.pkl'),'r'))
    
    # Pass the thermolin data into the model, get the corrections.
    # Return and print-to-screen the corrections to a5 and a6.
    # Make sure all issues of cos(dec) and sign are sorted out.
    x_predict = rf_x.predict(thermolin)
    y_predict = rf_y.predict(thermolin)

    assert(len(x_predict)==1)
    assert(len(y_predict)==1)
    x_predict = x_predict[0]
    y_predict = y_predict[0]

    # RK and JH believe the following sign convention is correct.
    add_to_col_x = +x_predict
    add_to_col_y = +y_predict

    # Add the corrections to the collimation terms.
    new_col_x = base_col_x + add_to_col_x
    new_col_y = base_col_y + add_to_col_y

    # Output the results to screen.
    print(' ')
    print('The size of the correction is (%+0.1f, %+0.1f) arcsec.'%(add_to_col_x*3600., add_to_col_y*3600.))
    print('The new collimate_fixed command is given below.  It goes into the ')
    print('online pointing config file, $ACTIVE_GCP_PATH/gcp/control/conf/spt/pointing.init.')
    print(50*'#')
    print('collimate_fixed eht, %0.4f, %0.4f'%(new_col_x, new_col_y))
    print(50*'#')
    print(' ')




def read_thermolin_data_from_archive(start_time, stop_time):
    '''
    This function fetches the (median) thermolin data between a 
    START_TIME and and STOP_TIME.  The times need to be datetime 
    objects.

    created by RK, 10 Jan 2015.
    '''

    from sptpol_software.data import arcfile
    registers = {'thermo':'antenna0.scu.temp', 
                 'lin':'antenna0.tracker.linear_sensor_avg~'}
    # NB: the linear sensor data is a "fast" register, while the thermometry 
    # data is not.  However, we're not going to keep track of any timestamps, 
    # because for the purposes of EHT pointing, all we're interested in 
    # is the median values of these registers during some short period of time.
    # In fact, this function is just going to return the median values.
    print('Getting thermolin data for ')
    print(start_time)
    print('to')
    print(stop_time)
    x = arcfile.readRegisters(registers, start_time, stop_time)
    med_thermo = np.median(x.thermo, axis=1)
    med_lin = np.median(x.lin, axis=1)
    # Now combine the median thermometry and linear sensor data into a single array.
    thermolin = np.hstack([med_thermo, med_lin])
    # There are 60 thermometry registers and 4 linear sensor registers.
    # The output numpy array will have 64 elements.
    return thermolin


def read_thermolin_data_from_archive_3G(start_time, stop_time, arcdir='/spt_data/arc'):
    pipe = core.G3Pipeline()
    pipe.Add(std_processing.ARCTimerangeReader, start_time=start_time, stop_time=stop_time, basedir=arcdir)
    pipe.Add(gcp.ARCExtract)
    
    d = {}
    d['linsens_avg_l1'] = []
    d['linsens_avg_l2'] = []
    d['linsens_avg_r1'] = []
    d['linsens_avg_r2'] = []
    d['scu_temp'] = []

    def fillPointingStructure(fr):
        if fr.type != core.G3FrameType.GcpSlow:
            return
        d['linsens_avg_l1'].append(fr['TrackerPointing'].linsens_avg_l1)
        d['linsens_avg_l2'].append(fr['TrackerPointing'].linsens_avg_l2)
        d['linsens_avg_r1'].append(fr['TrackerPointing'].linsens_avg_r1)
        d['linsens_avg_r2'].append(fr['TrackerPointing'].linsens_avg_r2)
        d['scu_temp'].append(fr['TrackerPointing'].scu_temp/core.G3Units.K+273.15)
    
    pipe.Add(fillPointingStructure)
    pipe.Run()

    d['linsens_avg'] = np.median(np.array([np.array(d['linsens_avg_l1']).flatten(), 
                                        np.array(d['linsens_avg_l2']).flatten(), 
                                        np.array(d['linsens_avg_r1']).flatten(), 
                                        np.array(d['linsens_avg_r2']).flatten()]),axis=1)
    d['scu_temp'] = np.median(np.array(d['scu_temp']),axis=0)

    thermolin = np.hstack([d['scu_temp'],d['linsens_avg']])

    return thermolin





def fit_thermolin_model(force_contiguous=False, do_save=False):
    '''
    This function is only meant to run "once".  Ryan K 
    is fitting pointing offsets measured by Jason H.
    They're doing this for EHT pointing, early Jan 2015.

    Fit a model to the pointing residuals using the 
    telescope thermometry and linear sensor data, aka 
    the "thermolin" data.

    There is a model for each direction (horizontal and 
    vertical), and we save each model to a pickle file for
    future use.
    '''

    # Load the thermolin data.  This includes
    # the "features" (FEAT) and the targets, the 
    # X and Y offsets we're trying to predict.
    data = load_residual_data(sources='hii')
    feat = data['thermolin']

    data_ptsrc = load_residual_data(sources='point')
    feat_ptsrc = data_ptsrc['thermolin']

    if do_save:
        data_all = load_residual_data(sources='all')
        feat_all = data_all['thermolin']

    name = {'xoff':'X', 'yoff':'Y'}
    for this_target in ['xoff','yoff']:
        target = data[this_target]
        target_ptsrc = data_ptsrc[this_target]

        if force_contiguous:
            nsamples = len(target)

            if False:
                # train of first 2/3, test on last 1/3
                ind = int(0.66*nsamples)
                feat_train = feat[0:ind, :]
                target_train = target[0:ind]
                feat_test = feat[ind:, :]
                target_test = target[ind:]
            else:
                # train on last 2/3, test on first 1/3
                ind = int(0.33*nsamples)
                feat_train = feat[ind:, :]
                target_train = target[ind:]
                feat_test = feat[0:ind, :]
                target_test = target[0:ind]           

        else:
            #from sklearn.cross_validation import train_test_split
            #feat_train, feat_test, target_train, target_test = train_test_split(feat, target, test_size=0.33)
            feat_train, feat_test, target_train, target_test = chunked_train_test_split(feat, target, test_size=0.33)

        # Fit a random forest regressor to the training data.
        from sklearn.ensemble import RandomForestRegressor
        rf = RandomForestRegressor(n_estimators=10, n_jobs=5)
        rf.fit(feat_train, target_train)

        # Make predictions for the redshifts of the test data.
        target_predict = rf.predict(feat_test)

        # Calculate the residual error
        err = (target_predict-target_test)
        print('HII   %s RMS: before %0.1f, after %0.1f, %0.2fX better'%(name[this_target], np.std(target)*3600., np.std(err)*3600., np.std(target)/np.std(err)))

        # And make predictions for the lower-elevation point sources.
        #print target_ptsrc
        target_predict_ptsrc = rf.predict(feat_ptsrc)
        err_ptsrc = (target_predict_ptsrc-target_ptsrc)
        print('PTSRC %s RMS: before %0.1f, after %0.1f, %0.2fX'%(name[this_target], np.std(target_ptsrc)*3600., np.std(err_ptsrc)*3600., np.std(target_ptsrc)/np.std(err_ptsrc)))

        # Now train using all data, and save the final model
        # for future use.
        if do_save:
            target_all = data_all[this_target]
            rf_all = RandomForestRegressor(n_estimators=10, n_jobs=5)
            rf_all.fit(feat_all, target_all)
            pickle.dump(rf_all, open('/home/rkeisler/eht_pointing/thermolin_model_%s.pkl'%(this_target),'w'))


def load_residual_data(remove_hii_median=True, 
                       sources=['rcw38','mat5a','0537-441','0521-365','2255-282'], 
                       years=['2013','2014'], 
                       max_ok_arcsec=100.):
    '''
    This function is only meant to run "once".  Ryan K 
    is fitting pointing offsets measured by Jason H.
    They're doing this for EHT pointing, early Jan 2015.

    Load the residual pointing data.
    '''
    residual_data_dir = '/data/jhenning/EHT_pointing/'
    from glob import glob
    if sources=='all': sources=['rcw38','mat5a','0537-441','0521-365','2255-282']
    if sources=='hii': sources=['rcw38','mat5a']
    if sources=='point': sources=['0537-441','0521-365','2255-282']
    output = {thing:[] for thing in ['xoff','yoff','thermolin']}
    el = {'rcw38':47.51, 'mat5a':61.33, '0537-441':44.085816, '0521-365':36.458569, '2255-282':27.97257}
    a4_fixed = -0.01358062
    a6_fixed = -0.21985015
    for source in sources:
        this_el = el[source]
        cos_el = np.cos(this_el*np.pi/180.)
        sin_el = np.sin(this_el*np.pi/180.)
        tan_el = sin_el/cos_el
        for year in years:
            filenames = sorted(glob('%s/%s/%s/*.pkl'%(residual_data_dir, source, year)))
            if len(filenames)==0: continue
            these_x = []
            these_y = []
            these_thermolin = []
            for filename in filenames:
                tmp = pickle.load(open(filename,'r'))
                # If this is the focus quasar 0537-441, only use observations taken during the 
                # "pointing" run in Nov/Dec 2014, as the other observations could be at a variety of 
                # focus positions.
                if source=='0537-441':
                    if (('2014112' not in filename) & ('2014120' not in filename)):
                        continue

                this_x = tmp['az_offset_divCosEl']
                this_y = tmp['el_offset_proj0']

                # Careful here.  Jason didn't apply the best-fit-global a4 and a6 when 
                # he made the maps that made these offsets, so I need to apply them here 
                # such that the thermolin correction has zero mean.
                this_x += a4_fixed*tan_el
                this_y -= a6_fixed

                # at this point the offsets are in "az/el" or "ra/dec" units.
                # let's escale them by some function of el, to see what 
                # form the thermolin correction looks most like. 
                # e.g. el tilt (a4) or cross-el collimation (a5)?
                #xsf = 1./tan_el # a4-like
                xsf = cos_el # a5-like
                # In fact, an el tilt (a4) performs slightly better than 
                # a cross-el collimation term (a5), but it's a small difference, 
                # and if we go with a5, then folks only have to change the 
                # online collimation (one line of pointing.init rather than 
                # two lines).
                ysf = 1.
                this_x *= xsf
                this_y *= ysf

                # go ahead and remove really crazy points.
                if np.sqrt( (this_x/xsf)**2. + (this_y/ysf)**2. ) > 2000: continue
                thermo = tmp['scu.temp']
                lin = tmp['tracker.linear_sensor_avg']

                # just use the thermometry and linear sensor data.
                this_thermolin = np.hstack([thermo,lin])

                '''
                # or, also include mbient temperature+pressure and wind.
                wind = np.array([tmp['observation.wind_speed_avg'], tmp['observation.wind_dir_avg']])
                ambient = np.array([tmp['observation.temp_avg'], tmp['observation.pressure_avg']])
                this_thermolin = np.hstack([thermo,lin,ambient,wind])
                '''

                these_x.append(this_x)
                these_y.append(this_y)
                these_thermolin.append(this_thermolin)

            if len(these_x)==0: continue

            these_x = np.array(these_x)
            these_y = np.array(these_y)
            these_thermolin = np.array(these_thermolin)

            # if desired, remove the median offsets per HII source per year.
            if ((remove_hii_median) & ((source=='rcw38') | (source=='mat5a'))):
                these_x -= np.median(these_x)
                these_y -= np.median(these_y)
            output['xoff'].append(these_x)
            output['yoff'].append(these_y)
            output['thermolin'].append(these_thermolin)

    output['xoff'] = np.hstack(output['xoff'])
    output['yoff'] = np.hstack(output['yoff'])
    output['thermolin'] = np.vstack(output['thermolin'])

    # remove outlying points.
    r = np.sqrt( (output['xoff']/xsf)**2. + (output['yoff']/ysf)**2. )
    whok = np.where(r<(max_ok_arcsec/3600.))[0]
    output['xoff'] = output['xoff'][whok]
    output['yoff'] = output['yoff'][whok]
    output['thermolin'] = output['thermolin'][whok,:]

    return output


def chunked_train_test_split(X, y, test_size=0.33, nchunks=100):
    '''
    If I have a set of time-ordered feature data X and label/target y, 
    I often want to split them up into two non-overlapping sets: 
    training and testing (or validation) sets, so get a sense of the
    accuracy of my classifier/regressor.  However, if the features and
    labels are correlated between neighboring time samples, you 
    will infer an artifically high accuracy.  This function gets around
    that by dividing the data into some numbers of contiguous chunks, 
    (each chunk contains multiple time samples), and then randomly 
    assigning the chunks to the training and testing sets.  This 
    mitigates the problem described above.
    '''
    nsamples = len(y)
    nsamples_per_chunk = int(np.floor(1.*nsamples/nchunks))
    ind_chunk = np.arange(nchunks)
    np.random.shuffle(ind_chunk)
    test_train_border = int(np.floor((1.-test_size)*nchunks))
    ind_chunk_train = ind_chunk[0:test_train_border]
    ind_chunk_test = ind_chunk[test_train_border:]

    # make training set
    X_train = []
    y_train = []
    for ind in ind_chunk_train:
        istart = ind*nsamples_per_chunk
        istop = (ind+1)*nsamples_per_chunk
        X_train.append(X[istart:istop, :])
        y_train.append(y[istart:istop])
    X_train = np.vstack(X_train)
    y_train = np.concatenate(y_train)

    # make testing set
    X_test = []
    y_test = []
    for ind in ind_chunk_test:
        istart = ind*nsamples_per_chunk
        istop = (ind+1)*nsamples_per_chunk
        X_test.append(X[istart:istop, :])
        y_test.append(y[istart:istop])
    X_test = np.vstack(X_test)
    y_test = np.concatenate(y_test)

    return X_train, X_test, y_train, y_test
