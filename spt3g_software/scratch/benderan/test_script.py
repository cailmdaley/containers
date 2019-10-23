from spt3g import core
frames=list(core.G3File('/spt/data/bolodata/downsampled/calibrator/3009897/0000.g3'))
ts=frames[-1]['RawTimestreams_I']['W136/2017.W136.5.4.7.912']
plt.plot(ts)



sync_sig=frames[-1]['CalibratorOn']

tmp=frames[2]['RawTimestreams_I']['W136/2017.W136.5.4.7.912']
sample_rate=tmp.sample_rate/core.G3Units.Hz   #this is in Hz now

tsec=np.arange(len(ts))/sample_rate



nbins=len(ts)
freq=np.array(range(0,(nbins-1)/2+2))*(sample_rate/nbins)
window=np.hanning(nbins)
fft_tmp=np.fft.fft(ts*window )
psd_tmp=np.zeros(len(freq))
psd_tmp[0]=np.real(fft_tmp[0]*np.conj(fft_tmp[0]))
psd_tmp[nbins/2]=np.real(fft_tmp[nbins/2]*np.conj(fft_tmp[nbins/2]))
psd_tmp[1:nbins/2]=2*np.real(fft_tmp[1:nbins/2]*np.conj(fft_tmp[1:nbins/2]))




