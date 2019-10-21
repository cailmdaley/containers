import numpy as np


## use distributions from data (eyeballed)
## data from .25' maps
T_MEAN  = {'90':140.,'150':150.,'220':180.}
T_STD   = {'90':8.,'150':13.,'220':20.}
NOISE_MEAN= {'90':3.4,'150':2.6,'220':10.6}
NOISE_STD = {'90':.7,'150':.4,'220':1.5}
BEAM_WIDTH = {'90':1.6,'150':1.3,'220':1.0}




def make_sim_CenA_TQU_template_maps():
    '''
    This function makes per-band T,Q,U maps (2D array) from 3 2D gaussians simulating CenA.
    I use polarization maps of data to estimate polarization angles, and amplitudes of lobes/nucleus.
    
    Arguments
    ---------
    
    band : str
        detector band
    polang : float
        Polarized emission angle (deg)
    polfrac : float
        Polarized emission fraction (in range [0.,1.])
        
    params for gaussian map maker : see make_2D_gaussian_Tmap
        defaults will use the BEAM_WIDTH[band] global variable
        and a resolution of .25'
        
    islobe : bool
        if True, reduces amplitude by 5x sin
        
    Returns
    -------
        
    T,Q,U : 2D numpy arrays
    '''
    ## define lobe/nucleus factors for amplitude, or polarization.
    upperlobefactor = {'90':0.35,'150':0.20,'220':.1}
    upperpolfrac    = {'90':0.25,'150':0.25,'220':.2}
    upperlobepolang = {'90':45.0,'150':45.0,'220':45.} 

    lowerlobefactor = {'90':0.2,'150':0.1,'220':.07}
    lowerpolfrac    = {'90':0.1,'150':.15,'220':.15}
    lowerlobepolang = {'90':90.,'150':90.,'220':90.}

    nuclearpolfrac = {'90':.03,'150':.06,'220':.07}#{'90':0.0,'150':0.0,'220':0.0} #
    nuclearpolang  = {'90':15.,'150':25.,'220':30.}
    
    bands = ['90','150','220']
    T = {}
    Q = {}
    U = {}
    for band in bands:
        Tnuc,Qnuc,Unuc = make_sim_TQU_maps(band,sigma=1.25,polfrac=nuclearpolfrac[band],polang=nuclearpolang[band])
        Tupperlobe,Qupperlobe,Uupperlobe = make_sim_TQU_maps(band,polfrac=upperpolfrac[band],polang=upperlobepolang[band],
                                                             sigma=2.6,mu=[-4,-4],ampfactor=upperlobefactor[band])
        Tlowerlobe,Qlowerlobe,Ulowerlobe = make_sim_TQU_maps(band,polfrac=lowerpolfrac[band],polang=lowerlobepolang[band],
                                                             sigma=2.,mu=[ 4,4 ],ampfactor=lowerlobefactor[band])

        T[band] = Tnuc+Tupperlobe+Tlowerlobe
        Q[band] = Qnuc+Qupperlobe+Qlowerlobe
        U[band] = Unuc+Uupperlobe+Ulowerlobe
    
    return T,Q,U
    
    
    
    
def make_2D_gaussian_Tmap(mapsize=120,sigma=1.0,mu=[0.,0.],res=.25):
    '''
    This function makes a 2D array with a 2D gaussian at the center by default.
    Width of the gaussian is symmetric in x,y and sigma is in units of arcmin.
    
    Arguments
    ---------
    
    mapsize : int
        side length of map (# of pixels)
    sigma : float
        width of gaussian (arcmin)
    mu : 2 element list
        mean of gaussian [x,y] (in arcmin)
    res : float
        map resolution (arcmin) -- converts the sigma to pixels
        
    Returns
    -------
    
    g : 2D numpy array
        2D array of size mapsize with a gaussian distribution defined above.
    
    '''
    sigma /= res
    mux = mu[0]/res
    muy = mu[1]/res
    x, y = np.meshgrid(np.arange(-mapsize//2,mapsize//2,1),np.arange(-mapsize//2,mapsize//2,1))
    d = np.sqrt((x-mux)**2+(y-muy)**2)
    g = np.exp(-( (d)**2 / ( 2.0 * sigma**2 ) ) )
    g[g<.001] = 0.
    return g

def make_sim_TQU_maps(band,polang=45,polfrac=.15,mapsize=120,sigma=-1,mu=[0.,0.],res=.25,ampfactor=1.):
    '''
    This function makes a T,Q,U map (2D array) from a 2D gaussian with T_MEAN amplitude.
    Q,U calculated with angle polang and polarized intensity fraction polfrac.
   
    Arguments
    ---------
    
    band : str
        detector band
    polang : float
        Polarized emission angle (deg)
    polfrac : float
        Polarized emission fraction (in range [0.,1.])
        
    params for gaussian map maker : see make_2D_gaussian_Tmap
        defaults will use the BEAM_WIDTH[band] global variable
        and a resolution of .25'
        
    islobe : bool
        if True, reduces amplitude by 5x sin
        
    Returns
    -------
        
    T,Q,U : 2D numpy arrays
        
    
    '''
    G      = 1.
    p      = 1.
    theta  = np.deg2rad(polang)
    ## Two equations for Q,U unknown
    ## eq(1): t = G * [(2-p)*T + p * cos(2Theta) * Q + p * sin(2Theta) * U]
    ## eq(2): np.sqrt(Q**2+U**2) = polfrac*T 
    
    if sigma == -1:
        T = make_2D_gaussian_Tmap(mapsize=mapsize,sigma=BEAM_WIDTH[band]/2.3,res=res,mu=mu)*T_MEAN[band]*ampfactor
    else:
        T = make_2D_gaussian_Tmap(mapsize=mapsize,sigma=sigma*BEAM_WIDTH[band]/2.3,res=res,mu=mu)*T_MEAN[band]*ampfactor
    
    
    ## pick positive root
    Q = np.cos(2*theta)*T*polfrac
    U = np.sin(2*theta)*T*polfrac
    
    return T,Q,U
    
    
    
def make_sim_bolo_map(T,Q,U,band,polfrac=.15,nomang=45,addnoise = False,randomize_noise=True,numobs=1.):
    '''
    This function takes simulated T,Q,U maps (2D arrays) and makes a simulated bolo map
    for detector with polarization angle nomang.
    
    Arguments
    ---------
    T,Q,U : 2D np.array
        simulated TQU maps
    band : str
        detector band
    nomang : float
        Polarized emission angle (deg)
    polfrac : float
        Polarized emission fraction (in range [0.,1.])
        
    params for gaussian map maker : see make_2D_gaussian_Tmap
        defaults will use the BEAM_WIDTH[band] global variable
        and a resolution of .25'
        
    Returns
    -------
        
    t : 2d-array
        simulated bolometer map
        
    
    '''
    G      = 1.
    p      = 1.
    theta  = np.deg2rad(nomang)
    ## eq(1): t = G * [(2-p)*T + p * cos(2Theta) * Q + p * sin(2Theta) * U]
    t = G * ((2-p)*T+ p * np.cos(2*theta) * Q + p * np.sin(2*theta) * U) #*(1-polfrac) 
    ## add noise 
    if addnoise:
        if randomize_noise:
            noisestd = np.random.normal(loc=NOISE_MEAN[band],scale=NOISE_STD[band])
        else:
            noisestd = NOISE_MEAN[band] 
        t += np.random.normal(loc=0.0,scale=noisestd/np.sqrt(numobs),size=np.shape(t))
    return t 
    
    
    
def make_sim_mask_fromT_map(T,cutamp = 1.):
    mask = np.zeros(np.shape(T))
    mask[np.where(T>cutamp)]+=1.
    return mask
    