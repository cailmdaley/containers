from spt3g import core
from spt3g.core import G3Units as U
from spt3g.util.genericutils import smooth, deriv, float_mode
import numpy as np

class AnalyzeElnod(object):

    def __init__(self,
                 bolos_to_use=None,
                 min_good_time=5 * U.s,
                 elrate_min=0.05 * U.deg/U.s,
                 max_difference_from_mode=0.1,
                 do_iq_rotation=True,
                 boresight_el_key='RawBoresightEl',
                 i_data_key='RawTimestreams_I',
                 q_data_key='RawTimestreams_Q'):
        self._clear()
        self.bolos_to_use = bolos_to_use
        self.min_good_time = min_good_time
        self.elrate_min = elrate_min
        self.max_difference_from_mode = max_difference_from_mode
        self.do_iq_rotation = do_iq_rotation
        self.boresight_el_key = boresight_el_key
        self.i_data_key = i_data_key
        self.q_data_key = q_data_key

    def _clear(self):
        self.bs_el = []
        self.i_data = []
        self.q_data = []

    def __call__(self, frame):

        # endprocessing signaled or elnod observation complete
        # and we have data to analyze
        if (len(self.bs_el) > 0 and 
            (frame.type == core.G3FrameType.EndProcessing
             or (frame.type == core.G3FrameType.Scan and 
                 'elnod' not in frame['GCPFeatureBits']))):

            # combine elnod frames into contiguous timestreams
            bs_el = core.G3Timestream.concatenate(self.bs_el)
            i_data = core.G3TimestreamMap.concatenate(self.i_data)
            if len(self.q_data) > 0:
                q_data = core.G3TimestreamMap.concatenate(self.q_data)
            else:
                q_data = None

            # do the thing
            cal_frame = analyze_elnod_data(
                bs_el, i_data, q_data,
                bolos_to_use=self.bolos_to_use,
                min_good_time=self.min_good_time,
                elrate_min=self.elrate_min,
                max_difference_from_mode=self.max_difference_from_mode,
                do_iq_rotation=self.do_iq_rotation)

            # clear and return
            self._clear()
            return [cal_frame, frame]

        # wrong type of frame
        if frame.type != core.G3FrameType.Scan:
            return

        # have not yet reached a turnaround
        if frame.get('Turnaround', False) and len(self.bs_el) == 0:
            return

        # we're in an elnod observation and after the first turnaround
        # so store this frame's data in the buffer
        if 'elnod' in frame['GCPFeatureBits']:
            self.bs_el.append(frame[self.boresight_el_key])
            self.i_data.append(frame[self.i_data_key])
            if self.q_data_key in frame:
                self.q_data.append(frame[self.q_data_key])

class RotateIQ(object):

    def __init__(self,
                 i_vector_key='ElnodEigenvalueDominantVectorI',
                 q_vector_key='ElnodEigenvalueDominantVectorQ',
                 i_data_key='RawTimestreams_I',
                 q_data_key='RawTimestreams_Q',
                 i_rotate_key='RawTimestreams',
                 q_rotate_key=None):
        self.iq_vec = None
        self.i_vector_key = i_vector_key
        self.q_vector_key = q_vector_key
        self.i_data_key = i_data_key
        self.q_data_key = q_data_key
        self.i_rotate_key = i_rotate_key
        self.q_rotate_key = q_rotate_key

    def __call__(self, frame):

        if frame.type == core.G3FrameType.Calibration:
            if self.i_vector_key in frame:
                self.iq_vec = [frame[self.i_vector_key], frame[self.q_vector_key]]

        elif frame.type == core.G3FrameType.Scan:
            if self.iq_vec is None:
                raise ValueError('Missing IQ phase calibration data')

            only_i = self.q_rotate_key is not None

            ret = rotate_iq(frame[self.i_data_key], frame[self.q_data_key], self.iq_vec,
                            only_i=only_i, inplace=False)
            if only_i:
                frame[self.i_rotate_key] = ret
            else:
                frame[self.i_rotate_key] = ret[0]
                frame[self.q_rotate_key] = ret[1]

def analyze_elnod_data( bs_el, 
                        i_ts_map, q_ts_map ,
                        bolos_to_use,
                        min_good_time, 
                        elrate_min, 
                        max_difference_from_mode,
                        do_iq_rotation):
    #check that data is good
    if len(i_ts_map.keys()) == 0:
        core.log_fatal("No data provided to analyze_elnod_data")

    #find the constant velocity portion of the elnod
    n_samples = len(bs_el)
    el_time_vals = np.linspace(0, n_samples / bs_el.sample_rate, n_samples, endpoint=False)
    el_velocity = smooth(deriv(el_time_vals, bs_el), 10, 'flat')
    vel_mode = float_mode(abs(el_velocity), n_bins = n_samples/50+1)
    vel_good = abs( abs(el_velocity) - vel_mode ) / vel_mode < max_difference_from_mode
    scanning_good = vel_good & (abs(el_velocity) > elrate_min)
    const_el_indices = np.where(scanning_good)[0]
    
    if len(const_el_indices) / bs_el.sample_rate < min_good_time:
        core.log_fatal("A long enough slice of constant el velocity was not found")

    '''
     We want a double regression on both time and airmass to account for linear drifts with time (like stage temp) and true airmass response.
     This is why a full elnod goes up and then down by the same amount -- to break this degeneracy.
     Correlation with time gives us drift.
     Correlation with airmass gives the sope, the elnod_response.
            
     With means subtracted from everything, the model is:
        this_data  =  drift * temp_time +  slope * airmass
     As a matrix problem, this is written as coefficients dotted into a template (again with means subtracted):
        this_data = [slope, drift]^T . [temp_time, airmass]
                  =          atemp^T . templates
     Least-squares solving (and defining X=templates),
         atemp = (X . X^T)^-1 . (X . data^T)
    '''

    el_data = bs_el[ const_el_indices ]
    airmass =  1.0 / np.sin( el_data / U.rad )
    temp_time = const_el_indices / bs_el.sample_rate
    templates = np.matrix( np.zeros( (2, len(el_data)) ) )
    templates[0] = temp_time - np.mean(temp_time)  # First  dependent variable is (mean subtracted) time
    templates[1] = airmass   - np.mean(airmass)    # Second dependent variable is (mean subtracted) airmass
    ctemp = np.linalg.inv(templates*templates.transpose()) # 2x2 matrix proportional to the inverse covariance matrix

    #initialize our output calibration frame
    frame = core.G3Frame(core.G3FrameType.Calibration)
    map_data_names = map(lambda s: 'Elnod'+s, 
                         ['EigenvalueRatio', 'EigenvalueDominantVectorI',
                          'EigenvalueDominantVectorQ',
                          'Slopes', 'Drifts', 'RSquared', 'Variances', 'SigmaSlopes','SNSlopes'] )
    for s in map_data_names:
        frame[s] = core.G3MapDouble()

    #loop through the bolos for actual analysis time
    for k in i_ts_map.keys():
        if q_ts_map is not None and not k in q_ts_map:
            core.log_fatal(k + " is in the i data but not the q data in the elnod analysis")

        i_data = i_ts_map[k][ const_el_indices ] 
        #iq analysis
        if q_ts_map is not None:
            q_data = q_ts_map[k][ const_el_indices ] 
            dom_eig_vec, sub_eig_vec, eig_ratio = get_dom_cov_eigenvec(i_data, q_data)
            if do_iq_rotation:
                rotate_iq(i_data, q_data, dom_eig_vec, only_i=True, inplace=True)

        #slope fitting analysis
        i_data -= np.mean(i_data)
        # 2x1 matrix. Correlate the data.  These are our fit coefficients
        atemp = ctemp * (templates * np.matrix(i_data).transpose()) 

        slope = float(atemp[1])  # airmass coefficient
        drift = float(atemp[0])  # time coefficient
        model = np.array(atemp.transpose() * templates).flatten()
        residuals = np.asarray(i_data - model)
        sst = np.sum((i_data - np.mean(i_data))**2)
        sse = np.sum(residuals**2)

        with np.errstate(divide='ignore', invalid='ignore'):
            r_squared = 1 - sse / sst
        variance = np.std(residuals)**2.0 #np.sum((residuals - np.mean(residuals))**2)/len(el_data)
        sigma_slope = np.sqrt(ctemp[1,1] * variance)
        sn_slope = 0.
        if np.abs(slope) > 0. and np.abs(sigma_slope) > 0.:
            sn_slope = np.abs(slope/sigma_slope)

        if slope > 0:
            '''
             Added 2013-Dec by Jason G
             Elnod response is the *wrong* sign.  In ADC or Watts, the correct sign is actually negative for the DfMUX convention.
             This means that the "headless" dominant eigenvector had the wrong sign and we have to rotate I/Q than +-90 degrees.
             The DfMUX people think this should never happen, but it happens for about one bolometer per observation.
             Fix everything here.
             We will be able to detect that this was done later because the rotation angle whose cos and sin are iq_projection[0] and [1] will be outside of +-90 deg.
            '''
            slope *= -1.
            drift *= -1.
            if q_ts_map is not None:
                dom_eig_vec[:] *= -1.  # both components of the vector.

        #fill the data into our calibration frame
        if q_ts_map is not None:
            frame['ElnodEigenvalueRatio'][k] = eig_ratio
            frame['ElnodEigenvalueDominantVectorI'][k] = dom_eig_vec[0]
            frame['ElnodEigenvalueDominantVectorQ'][k] = dom_eig_vec[1]

        frame['ElnodSlopes'][k] = slope
        frame['ElnodDrifts'][k] = drift
        frame['ElnodRSquared'][k] = r_squared
        frame['ElnodVariances'][k] = variance
        frame['ElnodSigmaSlopes'][k] = sigma_slope
        frame['ElnodSNSlopes'][k] = sn_slope

    return frame

#######################################################################################################

def get_dom_cov_eigenvec(i,q):
    '''
    i: i phase data
    q: q phase data

    returns domEigVec, subEigVec, eigRatio

    domEigVec: the normalized eigenvector with the largest eigenvalue  (hopefully close to [1,0])
    subEigVec: the normalized eigenvector with the smallest eigenvalue (orthogonal and hopefully close to [0,1])
    eigRatio: the ratio of eigenvalues. <= 1.  Tends to be 10^-4 at EL~40 and 10^-3 at EL~65

    If the eigenvalue ratio is tiny you are happy.
    If it's close to 1 it probably means this failed.
    If the timestream is all zero, return:
        domEigVec: [1,0],  subEigVec: [0,1],  eigRatio: nan
    
    Note that the input i and q can have huge DC offsets.
    These are effectively removed when we take the covariance matrix.
    It's also ok to apply the IQ rotation to timestreams with DC offsets --
    because the rotation is linear, 
    a huge DC offset in I and Q also gets rotated, 
    but importantly, the relative angle rotates by the same amount.
    
    HISTORY
        Created by Nick Harrington
        2014-Aug-07 JG - Handle case of zero timestreams.
            We used to return the same [0,1] eigenvector twice (wrong) and nan for the eigRatio (maybe ok)
            Now return "identity matrix" eigenvectors and still nan for the eigRatio
        2014-Dec-20 JG - Move to elnod.py since it made no sense to have this in readout.
    '''
    cm = np.cov(np.asmatrix([i,q]))
    
    eigVals, eigVecs = np.linalg.eigh( cm )

    if eigVals[0] >= eigVals[1]:
        domEigVec, subEigVec = eigVecs
        eigRatio = eigVals[1]/eigVals[0] if eigVals[0] != 0 else 0
    else:
        subEigVec, domEigVec = eigVecs
        eigRatio = eigVals[0]/eigVals[1] if eigVals[1] != 0 else 0

    # We want to preserve the sign of the I timestream in the data. 
    # Make sure the domEigVec has a positive first element.
    if domEigVec[0] < 0:
        domEigVec *= -1
        subEigVec *= -1  # probably doesn't matter, since these are only defined up to a sign anyway.
    #angle = np.deg2rad( np.arctan2(domEigVec[1], domEigVec[0]) )  # Angle to rotate i and q by.
    return domEigVec, subEigVec, eigRatio


##########################################################################################

def rotate_iq(i_data, q_data, domEigVec, only_i=False, inplace=False):
    '''
    Given the two-component unit vector domEigVec=(x,y), 
    which forms a small-ish angle with respect to the x-axis,
    rotate i_data and q_data between each other 
    as if each (i,q) pair was an (x,y) vector.
    Optionally do this *in place* to save time and memory.
    
       newI =  c*I + s*Q
       newQ = -s*I + c*Q
    
    modifies/returns I if only_i is True, otherwise returns (I, Q).
    '''
    # Get data with signal maximized for later processing. 
    # We expect the ratio of eigenvalues to be small if this worked.
    if inplace:
        I = i_data
        Q = q_data
    else:
        I = i_data.copy()
        Q = q_data.copy()
    c = domEigVec[0]  # cos(theta_IQ) to go into rotation matrix
    s = domEigVec[1]  # sin(theta_IQ)
    #print 'iq_rotation:', np.sqrt(c**2+s**2), c, s  # the magnitude is always 1.0, up to rounding.
    if only_i:
        # Just do the rotation for I.  Don't bother with Q.
        I[:] = c * I + s * Q
        return I
    # Do rotation for both I and Q. This has the form of a rotation matrix.
    (I[:], Q[:]) = (c * I + s * Q, -s * I + c * Q)
    return (I, Q)
