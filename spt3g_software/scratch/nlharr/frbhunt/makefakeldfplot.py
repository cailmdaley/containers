import numpy as np
import pylab as pl

xs = np.arange(2,51,1)


rates = 100.0/(xs)**0.5 + 20

pys = []
nys = []
for r in rates:
    pys.append(np.random.poisson(r))
    nys.append(np.random.poisson(20))

pys = np.array(pys)
nys = np.array(nys)

pl.errorbar(xs, pys, yerr=pys**0.5, fmt='.')
pl.errorbar(xs+0.2, nys, yerr=nys**0.5, fmt='.', color='r')
pl.xlabel("TES Separation (Arbitrary Units)")
pl.ylabel("Expected Counts for Pairs At Specific Separation")
pl.title("Expected Lateral Distribution Function Structure")
pl.legend(('Positive Events', 'Negative Events'))
pl.show()
