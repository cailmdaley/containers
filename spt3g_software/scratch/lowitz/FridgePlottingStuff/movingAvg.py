#moving average with window 

import numpy as np

def movingAvg(data,window):

    data=np.convolve(data, np.ones((window,))/window, mode='same')

    return data
