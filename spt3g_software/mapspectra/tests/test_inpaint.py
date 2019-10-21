#!/usr/bin/env python

from spt3g import core, coordinateutils
from spt3g.mapspectra import inpaint_map_laplace
import numpy as np
from copy import copy

####################
# inpainting constant
####################

xdim = 101
ydim = 201

m = coordinateutils.FlatSkyMap(xdim, ydim, core.G3Units.arcmin)
m += 10
mask = m.Clone(False)

msize = 10
np.asarray(mask)[ ydim//2-msize :ydim//2 + msize, xdim//2-msize :xdim//2 + msize] = 1
m *= (1 - mask)

inpaint_map_laplace(mask, 10000, m)

assert(np.max(np.abs(m - 10)) < 1e-5)

####################
# inpainting gradient
####################

xdim = 100
ydim = xdim

# create our gradient map
mg = np.meshgrid(np.arange(ydim), np.arange(xdim))
m = coordinateutils.FlatSkyMap(mg[0] + mg[1], core.G3Units.arcmin)
morig = copy(m)

# cut a hole
msize = 10
np.asarray(m)[ ydim//2-msize :ydim//2 + msize, xdim//2-msize :xdim//2 + msize] = 0

# make our mask
mask = m.Clone(False)
np.asarray(mask)[ ydim//2-msize :ydim//2 + msize, xdim//2-msize :xdim//2 + msize] = 1

# inpaint
inpaint_map_laplace(mask, 10000, m)

assert(np.max(np.abs((morig - m))) < 1e-5)
