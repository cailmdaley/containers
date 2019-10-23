import os
import numpy as np
from spt3g import core, std_processing

'''
    Explanation of function
    
    Arguments
    ---------
    arg name : arg type
        explanation of arg
        
    Returns
    -------
    output name : output type
        explanation of output
'''



def load_bolometer_properties(cal_file):
    '''
    loads bolo properties to dictionaries keyed by bolo.
    
    Arguments
    ---------
    cal_file : str
        path to the calibration file where bolometer properties are stored.
        
    Returns
    -------
    bolo2band,pname,AB,nominal_angles,coupling: dictionaries
        dictionaries keyed by bolometer which contain the information for which they're named.
    unique_nominal_angles : list
        list of the unique nominal angles
    '''

    bolo2band = {}
    pname     = {}
    AB        = {}
    coupling  = {}
    nominal_angles = {}
    unique_nominal_angles = []
    
    offline_cal = cal_file
    for frame in core.G3File(offline_cal):
        for bolo in frame['BolometerProperties'].keys():
            try:
                bolo2band[bolo]=int(frame['BolometerProperties'][bolo].band/10.0)
                if bolo in pname:
                    if frame['BolometerProperties'][bolo].physical_name != pname[bolo]:
                        print("Not same pname!")
                        print(bolo,pname[bolo],frame['BolometerProperties'][bolo].physical_name,ob)
                   
                pname[bolo] = frame['BolometerProperties'][bolo].physical_name
                AB[bolo] = frame['BolometerProperties'][bolo].pixel_type
                nomang = int(frame['BolometerProperties'][bolo].pol_angle/core.G3Units.deg)
                coupling[bolo] = frame['BolometerProperties'][bolo].coupling
                
                ## only care about range 0-180, since pol shift of 180 is degenerate
                if nomang< 0:
                    nomang+=180
                if nomang > 180:
                    nomang-=180
                    
                
                nominal_angles[bolo] = nomang
                
                if nominal_angles[bolo] not in unique_nominal_angles:
                    unique_nominal_angles.append(nominal_angles[bolo])
            except ValueError:
                pass
    
    return bolo2band,pname,AB,nominal_angles,unique_nominal_angles,coupling


def get_cena_obsids(data_loc = '/spt/user/production/calibration/calframe/CenA-pixelraster/',season='any',years = [],obs2exclude=[]):
    '''
    Globs for the observations in the CenA calframe data folder (or the data_loc pointed to)
    
    Arguments
    ---------
    data_loc : str
        path to the folder where the calframes are stored
        
    years : list
        years of the observations you want to find obsIDs for
        Will glob 4*, 5* for 2018, and 6*,7* for 2019, as of 5/21/2019
    
    Returns
    -------
    obsids: list
        list of the obsids located in data_loc
        
    '''
    
    import glob
    if years == []:
        data_files =glob.glob(data_loc+'4*.g3') + glob.glob(data_loc+'5*.g3') + glob.glob(data_loc+'6*.g3')+ glob.glob(data_loc+'7*.g3')
    elif 2018 in years:
        data_files = glob.glob(data_loc+'4*.g3') + glob.glob(data_loc+'5*.g3')
    elif 2019 in years:
        data_files = glob.glob(data_loc+'7*.g3')+glob.glob(data_loc+'8*.g3') ## could add 6* in, but that's sun-up. glob.glob(data_loc+'6*.g3') + 
    obsids = []
    for file in data_files:
        obs = file.split('/')[-1].split('.')[0]
        if season !='any':
            if int(obs)>70000000:
                continue
        if obs not in obs2exclude:
            obsids.append(obs)
    return obsids 
       
    
    
def remove_outliers(data, m=2):
    '''
    Removes m-sigma outliers from the median of a set of data
    
    returns original dataset if length of new array not > 1
    
    Arguments
    ---------
    data : list or array of data
    
    Returns
    -------
    data: numpy array
        data with m-sigma outliers removed
    
    '''
    data = np.asarray(data)
    #print(data)                                                                                                                                 
    gooddata = np.asarray(abs(data - np.nanmedian(data)) < m * np.nanstd(data))
    #print(gooddata)                                                                                                                             
    if len(gooddata)>1:
        return data[gooddata]
    else:
        return data


def calc_median_and_1sigma(data,sigmaclip=3):
    '''
    Takes in an array and returns the median and standard deviation.
    Option to remove outliers; if sigmaclip = 0, no clipping is performed.
    
    Arugments
    ----------
    
    data : list or array 
        data to be averaged
    
    sigmaclip : int
        +/- n sigma outside which to remove outliers.
    
    Returns
    ---------
    med : float
        Median of data
    
    sigma : float
        1-sigma from median of data
    
    '''
    data = np.array(data)
    if sigmaclip>0:
        data = remove_outliers(data,m = sigmaclip)
    
    med   = np.nanmedian(data)
    sigma = np.nanstd(data)
    
    return med,sigma
        
        
def make_directory(path):
    '''
    Checks for the existance of a path, and creates a folder within that path
    
    Arguments
    ---------
    path : str
        path to directory
     
    Returns
    -------
    
    
    '''

    if not os.path.exists(path):
        os.makedirs(path)
    
    return


def get_human_readable_time_from_seconds(t):
    '''
    Converts a timestamp in seconds to a human readable string.
    
    Arguments
    ---------
    t : float or str
        time in seconds
    
    Returns
    -------
    str, a readable string : [h]hours:[m]minutes:[s]seconds
    
    '''
    t = int(t)
    h = int(t/3600)
    m = int((t-(h*3600))/60)
    s = int(t-(h*3600)-(m*60))
    return '[%i]hours:[%i]minutes:[%i]seconds'%(h,m,s)
