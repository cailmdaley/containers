import numpy as np
import pylab as py
import scipy.optimize as so
import sptpol_software.util.mpfit as mpfit
import sptpol_software.analysis.pointing.pointing_tools as pt
import gather_offsets as go

py.ion()


year = '2016'

snr_cut = 15.
offset_cut = 15.
hit_cut = 0 #5.
#offset_dir = '/data55/bleeml/SPTpol_500d_Nov16/pointing_checks/2016/'
offset_dir = '/data55/bleeml/SPTpol_500d_Nov16/pointing_checks/2016_try3/'
search_str = 'smap_weighted*txt'

#Start and stop files to gather offsets from
file_start = 2
file_stop = None

fit_for_boom_terms = True



#CORRECTION model
def model_daz(el, a4, a5, a7, az0):
    d_az = a4*np.tan(el*np.pi/180.) - a5/np.cos((el-a7)*np.pi/180.) - az0

    return d_az

#CORRECTION model
def model_del(el, a0, a1, a6):
    d_el = a0*np.sin(el*np.pi/180.) + a1*np.cos(el*np.pi/180.) - a6

    return d_el

def mpfitfun(el,y,err):
    def f(p, fjac=None):
        a0 = p[0]
        a1 = p[1]
        a6 = p[2]

        model = model_del(el, a0, a1, a6)
        
        return [0,(model - y)/err]

    return f

def mpfitfun_az(el,y,err):
    def f(p, fjac=None):
        a4 = p[0]
        a5 = p[1]
        a7 = p[2]
        az0 = p[3]

        model = model_daz(el, a4, a5, a7, az0)
        
        return [0,(model - y)/err]

    return f

#Expected positions of RCW38 and MAT5A read in from sources.py
rcw38_position = np.array([134.77, 47.51])
mat5a_position = np.array([167.88625, 61.36222222])

deltas, az, el, az_offset, el_offset, errors, offsets, jitter = go.gather_offsets(offset_dir, snr_cut=snr_cut, 
                                                                                  hit_cut=hit_cut, 
                                                                                  offset_cut=offset_cut,
                                                                                  search_str=search_str,
                                                                                  file_start=file_start,
                                                                                  file_stop=file_stop)

#az_offset *= np.cos(el*np.pi/180.)


#Get the daz model
pa = []
for j in range(4):
    pa.append({'value':1, 'fixed':0, 'mpside':2, 'relstep':0.01,
               'limited':[0,0], 'limits':[-1.0,1.0]})

a4_guess = -0.45
a5_guess = -0.143
a7_guess = 0.0
az0_guess = -0.304527 #Best-fit from optical measurements.

#a4
pa[0]['value'] = a4_guess
pa[0]['limited'] = [1,1]
pa[0]['limits'] = [-0.75, 0.75]

#a5
pa[1]['value'] = a5_guess
pa[1]['limited'] = [1,1]
pa[1]['limits'] = [-1.0, 1.0]
pa[1]['fixed'] = 0

#a7 #We are no longer using this.  Keep fixed at 0.0.
pa[2]['value'] = a7_guess
pa[2]['limited'] = [1,1]
pa[2]['limits'] = [-25, 25]
pa[2]['fixed'] = 1

#az0
pa[3]['value'] = az0_guess
pa[3]['limited'] = [1,1]
pa[3]['limits'] = [-1., 1.]
pa[3]['fixed'] = 1

result_az = mpfit.mpfit(mpfitfun_az(el=el, y=az_offset, err=errors), 
                        parinfo=pa,xtol=1e-10,maxiter=5000,quiet=True)

popt_az = [result_az.params[0], result_az.params[1], result_az.params[2], result_az.params[3]]

test_el = (np.arange(36000)*0.01)

pa = []
for j in range(3):
    pa.append({'value':1, 'fixed':0, 'mpside':2, 'relstep':0.01,
               'limited':[0,0], 'limits':[-1.0,1.0]})

a0_guess = 0.00549382
a1_guess = -0.01758231
a6_guess = 0.26

#a0
pa[0]['value'] = a0_guess
pa[0]['limited'] = [1,1]
pa[0]['limits'] = [-0.75, 0.75]
pa[0]['fixed'] = int(fit_for_boom_terms)

#a1
pa[1]['value'] = a1_guess
pa[1]['limited'] = [1,1]
pa[1]['limits'] = [-0.75, 0.75]
pa[1]['fixed'] = int(fit_for_boom_terms)


#a6
pa[2]['value'] = a6_guess
pa[2]['limited'] = [1,1]
pa[2]['limits'] = [-0.75, 0.75]
pa[2]['fixed'] = 0

result_el = mpfit.mpfit(mpfitfun(el=el, y=el_offset, err=errors), 
                        parinfo=pa,xtol=1e-10,maxiter=5000,quiet=True)

print('a0, a1, a6 = ', result_el.params)
print('a4, a5, a7, az0 = ', result_az.params[0], result_az.params[1], result_az.params[2], result_az.params[3])

print('\n')

test_yel = model_del(test_el, result_el.params[0], result_el.params[1], 
                              result_el.params[2])
test_yaz = model_daz(test_el, result_az.params[0], result_az.params[1], 
                              result_az.params[2],result_az.params[3])


az_residuals = az_offset - model_daz(el, result_az.params[0], result_az.params[1], 
                                         result_az.params[2], result_az.params[3])
el_residuals = el_offset - model_del(el, result_el.params[0], result_el.params[1], result_el.params[2])


print('HII position corrections: \\n')
print('Year: ', year)

#plot data.
py.figure(1)
py.clf()
py.plot(el, az_offset, 'ko')
py.plot(test_el, test_yaz, 'k-')
py.xlabel('elevation')
py.ylabel('Az offset')
py.title('Az offsets: az0 = ' + str(result_az.params[3]))
py.ylim((-0.15,0.15))
py.xlim((45, 70))

py.figure(2)
py.clf()
py.plot(el, el_offset, 'ko')
py.plot(test_el, test_yel, 'k-')
py.xlabel('elevation')
py.ylabel('el offset')
py.title('El offsets: az0 = ' + str(result_az.params[3]))
py.xlim((45, 70))
py.ylim((0.15,0.33))

py.figure(3)
py.clf()
py.plot(el, az_residuals*3600., 'ko')
py.xlabel('elevation')
py.ylabel('Az residuals [arcsec]')
py.title('Az Residuals: az0 = ' + str(result_az.params[3]))

py.figure(4)
py.clf()
py.plot(el, el_residuals*3600., 'ko')
py.xlabel('elevation')
py.ylabel('El residuals [arcsec]')
py.title('El Residuals: az0 = ' + str(result_az.params[3]))


