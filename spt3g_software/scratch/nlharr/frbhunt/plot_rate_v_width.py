import pylab as pl
import numpy as np
import scipy.interpolate


data = [[1, 58417.06113081337, 128418.94029178804],
        [5, 107714.15913149976, 239337.4107933324],
        [10, 197434.877492749, 506034.71097704576]]

cs = [(30/256., 191/256., 216/256.),
      (102/256., 74/256., 158/256.),
      (206/256., 26/256., 170./256)]
xs, one_sig, two_sig = zip(*data)

for i, ys in enumerate([one_sig, two_sig][::-1]):
    ys = np.array(ys)
    '''
    pl.fill_between(xs, ys*0, ys, facecolor=cs[1-i],
                    linewidth=0, alpha = 1.0)
    pl.plot(xs, ys, color = cs[1-i])
    '''
    xs = np.array(xs)
    pl.errorbar(xs+i*0.0, 0*xs, yerr = [0*xs, ys], fmt = '.', linewidth = 15, color=cs[1-i])

pl.xlim([0,11])
pl.ylim([0,7e5])
lgnd = pl.legend(['0.68 Confidence', '0.95 Confidence'], numpoints=1)
lgnd.legendHandles[0]._legmarker.set_markersize(15)
lgnd.legendHandles[1]._legmarker.set_markersize(15)

pl.xlabel("Input Signal FWHM (ms)")
pl.ylabel('Rate (Transient Events $sky^{-1} day^{-1})$ Above 10 Jy ms')
pl.title("Sky Rate Constraints as a Function of Signal Time FWHM")
pl.show()

