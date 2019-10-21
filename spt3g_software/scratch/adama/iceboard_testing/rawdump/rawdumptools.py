# rawdumptools.py
#
# Python script to process the data from the high-frequency FPGA "raw dump".
# We process the data file by file, breaking each file into smaller time blocks
# for which we compute the PSD. The PSD data and plots are saved and the
# original raw data files (which are very large) are deleted to save space.
# There is an option to keep raw data.
#
# Adam Anderson
# adama@fnal.gov

import sys, shutil
import numpy as np
import cPickle as pickle

def makePSD(ind_start_average, fft_length, n_fft_average):
    freq = scipy.fftpack.rfftfreq(fftlen, 1./20e6)
    psd = numpy.zeros(len(freq))
    j_fft = 0
    while j_fft < n_fft_average:
        ind_start_fft = ind_start_average + (j_fft*fft_length)
        fft_sample = scipy.fftpack.rfft(data[ind_start_fft:(ind_start_fft + fft_length)])
        psd += (numpy.conj(fft_sample)*fft_sample).real / n_fft_average

    return freq, psd


def process(infilename, outfilename, fft_length=65536, n_fft_average=300, delete=False):
        # load the data as a 'memmap' to reduce memory consumption
        data = numpy.memmap(infilename, np.int16, 'r')

        # compute the PSD of the data in chunks of length fft_length and average
        # N chunks together. repeat this process until we reach the end of the file
        
        psds = []
        ind_start_average = 0
        while ind_start_average < (len(data) - (n_fft_average*fft_length)):
            freq, psd = makePSD(ind_start_average, fft_length, n_fft_average)
            psds.append(psd)
            ind_start_average += (n_fft_average*fft_length)

        # save averaged PSDs as pickle files
        psd_data = {'freq': freq, 'PSDs': np.asarray(psds)}
        processed_file = open(outfilename, 'w')
        pickle.dump(psd_data)
        processed_file.close()

        # delete raw data file if requested
        os.remove(infilename)
        
