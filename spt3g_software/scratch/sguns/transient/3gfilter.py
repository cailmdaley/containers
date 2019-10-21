from spt3g import core, mapmaker,coordinateutils
import numpy as np
import tmplfilter
import sys,re,argparse


def gaussian(height, center_x, center_y, width_x, *args):
    width_x = float(width_x)
    width_y = float(width_x)
    rotation = 0
    if len(args) > 0:
        width_y = float(args[0])
    if len(args) > 1:
        rotation = float(args[1])
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2 + 2*rotation*(center_x-x)*(center_y-y))/2)/(2*np.pi*width_x*width_y)


def filterframe(frame):
    if frame.type != core.G3FrameType.Map:
        return
    if not ('Id' in frame and 'GHz' in frame['Id']):
        return

    band = re.search(r"[0-9]{2,3}GHz", frame['Id']).group(0)
    res_arcmin = frame['T'].res / core.G3Units.arcmin
    alpha_center = frame['T'].alpha_center
    delta_center = frame['T'].delta_center
    proj = frame['T'].proj
    # Conversion to mJy:
    # - Use SPT-SZ data release paper table 1 to get 396.3e6 Jy/K/sr
    # - Convert sr to arcmin^2 (8.46159496e-8)
    # - Convert Jy to mJy (10^3)
    # - Multiply by resolution^2 to get per-pixel density

    tomjy={'90GHz': 202.7e6*8.46159496e-8*1000*res_arcmin**2,#!!!: using delta function band approx for 90GHz
           '150GHz': 396.3e6*8.46159496e-8*1000*res_arcmin**2,
           '220GHz': 476.7e6*8.46159496e-8*1000*res_arcmin**2}

    filter = np.zeros((45,45)) 
    r = np.hypot(np.indices(filter.shape)[0] - filter.shape[0]/2, np.indices(filter.shape)[1] - filter.shape[1]/2)
    filter[r > 2./res_arcmin] = 1
    filter[r > 5./res_arcmin] = 0

    if 'Wpol' in frame:
        tuw, quw, uuw = mapmaker.remove_weight(frame['T'], frame['Q'], frame['U'], frame['Wpol'])
        del quw, uuw
        frame.pop('T',0); frame.pop('Q',0); frame.pop('U',0); 
        csigma = np.sqrt(1./frame['Wpol'].TT)
        frame.pop('Wpol',0)

    elif 'Wunpol' in frame:
        tuw = mapmaker.mapmakerutils.remove_weight_t(frame['T'], frame['Wunpol'])
        frame.pop('T',0);
        csigma = np.sqrt(1./frame['Wunpol'].TT)
        frame.pop('Wunpol',0)

    tuw /= core.G3Units.K
    csigma /= core.G3Units.K

    #print "Filtering..."
    filtered,filteredsigma = tmplfilter.tmplfilter(tuw, csigma, filter)
    #print "Done with Real Space hi ell filter."

    cm = (tuw - filtered)*tomjy[band]
    del tuw, filtered
    csigma = np.hypot(csigma, filteredsigma)*tomjy[band]
    del filteredsigma

    csigma *= np.std((cm/csigma)[np.isfinite(cm/csigma)])


    #print "Point Source filtering..."
    beamsigma = {'90GHz':1.54/2.35482/res_arcmin, '150GHz': 1.13/2.35482/res_arcmin, '220GHz':1.0/2.35482/res_arcmin}
    filterwidth = int(np.ceil(3*beamsigma[band]))
    psfilter = gaussian(1, filterwidth, filterwidth, beamsigma[band])(*np.indices((2*filterwidth+1,2*filterwidth+1)))
    psfiltered, psfilteredsigma = tmplfilter.tmplfilter(cm, csigma, psfilter)
    del cm, csigma
    #print "Done with Point Source filter."

    outframe = core.G3Frame(core.G3FrameType.Map)
    outframe['Tfilt'] = mapmaker.mapmakerutils.FlatSkyMap(psfiltered, res = res_arcmin*core.G3Units.arcmin, proj = proj, alpha_center = alpha_center, delta_center = delta_center)
    outframe['Sigmafilt'] = mapmaker.mapmakerutils.FlatSkyMap(psfilteredsigma, res = res_arcmin*core.G3Units.arcmin, proj = proj, alpha_center = alpha_center, delta_center = delta_center)
    outframe['Band'] = band
    return outframe

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', action ='store')
    parser.add_argument('-o','--output', action ='store')
    pargs = parser.parse_args()

    p = core.G3Pipeline()
    p.Add(core.G3Reader, filename = pargs.input)
    p.Add(filterframe)
    p.Add(lambda fr: 'Tfilt' in fr)
    p.Add(core.G3Writer, filename = pargs.output)
    p.Run(profile=True)
