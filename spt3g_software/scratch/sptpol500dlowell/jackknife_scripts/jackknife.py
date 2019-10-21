from spt3g import core, coordinateutils, mapmaker, mapspectra
from spt3g.mapspectra import map_analysis
import numpy, pickle, sys, copy
import glob
import argparse as ap
import numpy as np 
from random import shuffle
'''
Makes a jackknife from two arguement sets
'''

P = ap.ArgumentParser(description='Jackknife Bundles of data',
              formatter_class=ap.ArgumentDefaultsHelpFormatter)
P.add_argument('input_files', action='store', nargs='+', default=[], help='two pickle files with lists of the maps to go into the two groups')
P.add_argument('-a', '--apod', action='store', default='/spt/user/javva/lowell/good_apodization_mask.pkl', help='Apodization Mask')
P.add_argument('-eb', '--eb', action = 'store_true', default = False)
P.add_argument('-ron', '--read_out_nulls', action = 'store', default = None)
P.add_argument('-npl', '--plot', action = 'store_false', default = True)
P.add_argument('-o', '--output', action='store', default='output.pkl',
           help='Output filename')


args = P.parse_args()

group_1 = pickle.load(open(args.input_files[0], "rb"))
group_2 = pickle.load(open(args.input_files[1], "rb"))

if not args.eb:
        qu = True 

#loads in apodization mask -- note default is for sptpol 500d
apod = pickle.load(open(args.apod, "rb"))
numpy.asarray(apod)[numpy.asarray(apod) < 0] = 0

#binning of PS
c = numpy.arange(25, 5000, 25)
ell_bins = mapspectra.basicmaputils.get_reg_spaced_ell_bins(c)

def frame_from_uwmaps(t,q,u):
	f = core.G3Frame()
	f['T'] = t
	f['Q'] = q
	f['U'] = u
	return f

def explode_ps_dict(d):
	if qu:
		return d['TT'], d['QQ'], d['UU']
	else:
		return d['TT'], d['EE'], d['BB']
#stacks the two groups
g1 = []
print('These are the files in group 1')
for i in group_1:
        for f in core.G3File(i):
                print(i)
                g1.append((f['T'], f['Q'], f['U'], f['Wpol']))
g2 = []
print('These are the files in group 2')
for i in group_2:
	for f in core.G3File(i):
		print(i)
		g2.append((f['T'], f['Q'], f['U'], f['Wpol']))

print('Creating the null maps')
#makes null maps from pairs
nulls = []
for idx in range(len(g1)):
	print('Making null number' , idx )
	a = mapmaker.mapmakerutils.remove_weight(*g1[idx])
	b = mapmaker.mapmakerutils.remove_weight(*g2[idx])
	dt, dq, du = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
	nulls.append((dt, dq, du))

if args.read_out_nulls:
        print('Reading out null maps to %s'%args.read_out_nulls)
        for idx, i in enumerate(nulls):
                f = core.G3Frame(core.G3FrameType.Map)
                f['T'] = i[0]
                f['Q'] = i[1]
                f['U'] = i[2]
                core.G3Writer(args.read_out_nulls+'null_%s.g3'%idx)(f)
        nulls_info = {}
        nulls_info['group1'] = group_1
        nulls_info['group2'] = group_2
        with open(args.read_out_nulls+'nulls_defs.pkl', 'wb') as handle:
                pickle.dump(nulls_info, handle, protocol=pickle.HIGHEST_PROTOCOL)

noisecls = []
xcls = []
pairs = [(i, j) for i in range(len(nulls)) for j in range(len(nulls)) if i != j]

#calculates cross spectra for all the possible pairs of null maps
for i,j in pairs:
        print(' Calc xspectra from:', i,j)
        f1 = frame_from_uwmaps(*nulls[j])
        f2 = frame_from_uwmaps(*nulls[i])
        crosst, crossq, crossu = explode_ps_dict(map_analysis.calculateCls(f1,f2, apod_mask=numpy.asarray(apod), qu=qu, ell_bins=ell_bins))
        xcls.append((crosst, crossq, crossu))

#average crosspectra and error bar	
xcls = numpy.asarray(xcls)
meanpower = numpy.mean(xcls, axis=0)
t,q,u = meanpower
noisestdpower = numpy.std(xcls, axis=0)/numpy.sqrt(len(xcls))

output_file = {}
output_file['xcls'] = xcls
output_file['mean_power'] = {}
output_file['mean_power']['T'] = t
output_file['mean_power']['Q'] = q
output_file['mean_power']['U'] = u
output_file['noise_std_power'] = noisestdpower

with open(args.output, 'wb') as handle:
        pickle.dump(output_file, handle, protocol=pickle.HIGHEST_PROTOCOL)


if args.plot:
        pylab.figure()
        # XXX units incorrect
        #pylab.errorbar(c, t *c*(c+1), noisestdpower[0]*c*(c+1), label='T')
        pylab.errorbar(c, q / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1), noisestdpower[1] / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1), label='Q')
        pylab.errorbar(c, u / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1), noisestdpower[2] / (core.G3Units.arcmin * core.G3Units.uK)*c*(c+1),label='U')
        pylab.semilogx()
        pylab.legend()
        pylab.xlabel('Ell')
        pylab.ylabel('C_ell')
        pylab.grid()

