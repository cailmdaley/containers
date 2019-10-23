import numpy as np
import sys
sys.path.append('/home/axf295/2019/code/spt3g_software/scratch/axf295/polcal/python/')
import Circle_Functions as CF
from matplotlib import pylab as plt
import time
from scipy.stats import chi2 as scipy_chi2





def rotate_nomang_fits(nomangfits,rot_ang):
    ## Rotate all the fits (for all wafers and the coadd)
    ## uses the best fit rotation angle 
    wafers = nomangfits.keys()
    for wafer in wafers: 
        for ang in nomangfits[wafer]:
            for band in nomangfits[wafer][ang]:
                phiold = np.copy(nomangfits[wafer][ang][band]['rotPhi'])
                phinew = phiold +rot_ang
                if phinew < 0.:
                    phinew+=180.
                if phinew > 180.:
                    phinew -= 180.
                nomangfits[wafer][ang][band]['rotPhi'] = phinew
                nomangfits[wafer][ang][band]['res']  = CF.difference(phinew,ang,deg=True)
                
    return nomangfits
   
def get_optimal_nomang_rotation(nomangfits):
    ## minimize chi^2 by rotating nominal angles
    ## fit to 90GHz since those are the highest SN
    nomangs = sorted(nomangfits.keys())
    fits = []
    stds = []
    for ang in nomangs:
        fits.append(nomangfits[ang]['90']['rawPhi'][0])
        stds.append(nomangfits[ang]['90']['rawPhi'][1])
    var = np.var(np.array(fits)-nomangs)
    dof = len(nomangs)
    rotangle = np.linspace(0,180.,181)
    chi2 = []
    for i in range(len(rotangle)):
        angs = np.copy(nomangs)-rotangle[i]
        angs[angs>180.]-=180.
        angs[angs<0.]+=180.
        chi2.append(np.nansum((CF.difference(angs,fits,deg=True))**2/stds)/(dof-1))
    bestang = rotangle[np.where(chi2==min(chi2))]
    return bestang[0]

        
    

def optimize_map_alignment(img,background,band,xshift,yshift):
    ## Shift img map to find shift which minimizes the rms
    ## residual. img is shifted while background remains fixed.
    ## map1 and map2 need to have same shape
    xrange = np.arange(-1*xshift,xshift+1,1)
    yrange = np.arange(-1*yshift,yshift+1,1)
    ## assuming square image
    s = np.shape(img)[0]//2
    boundary  = s - np.max([xshift,yshift])
    
    res = np.zeros((len(xrange),len(yrange)))
    x = 0
    for i in xrange:
        y = 0
        for j in yrange:
            #print(i,j,x,y)
            #print(s-boundary+i,s+boundary+i,s-boundary+j,s+boundary+j,s-boundary,s+boundary)
            res[x][y] = np.nansum(abs(img[s-boundary+i:s+boundary+i,s-boundary+j:s+boundary+j]-background[s-boundary:s+boundary,s-boundary:s+boundary]))
            y+=1
        x+=1
    #plt.figure();plt.imshow(res);plt.colorbar();plt.savefig('/home/axf295/2019/tmpPlots/%sGHz_mapalignment.png'%band);plt.close()
    minvals = np.where(res == np.amin(abs(res)))
    if len(minvals[0])==0:
        minvals = np.where(res == -np.amin(abs(res)))
    fp = [xrange[minvals[0][0]],yrange[minvals[1][0]]]
    img = np.pad(img[s-boundary+fp[0]:s+boundary+fp[0],s-boundary+fp[1]:s+boundary+fp[1]],
                 (xshift,yshift),'constant')
    print(fp) 
    return img


def normalize_TQU_template_maps(tqu_template_maps, perband_T_coadds,savedir = './', radius = 1.,res = .25,source='CenA',obs='', debug = False):
    '''
    Provides a normalization to TQU template maps so that the coadd T map
    amplitude is equal to the TQU template amplitude.
    
    Input maps are assumed to be centered on the source. A small 1' radius around the
    center is used by default to gain-match the maps.
    
    Arguments
    ---------
    tqu_template_maps : dict
        Template maps used to do the fits; dictionary is keyed by band and contains CenAMap classes.
        
    
    perband_T_coadds : dict
        coaddmaps made from coadding single bolometer maps. dictionary of 2D numpy arrays
        Expected that this should be equal to T of the TQU template... but doesn't seem to be that way
        hence this function.
        
    radius : float
        radius, in arcmin, to use as the gain-matching region.
        
    res : float
        map resolution in arcmin.
        
    Returns
    -------
    normalized tqu_template_maps : dict
        TQU template maps divided by the factor which minimizes the difference between the coadded single bolo T
        map and the T template map.
    '''
    from scipy.interpolate import griddata
    
    dr = int(radius/res)
    
    for band in perband_T_coadds:
        if band in tqu_template_maps:
            t = perband_T_coadds[band]
            
            T = tqu_template_maps[band].maps['T']
            
            tcenter = np.shape(t)[0]//2
            Tcenter = np.shape(T)[0]//2
                       
            fluxt = np.sum(t[tcenter-dr:tcenter+dr,tcenter-dr:tcenter+dr])
            fluxT = np.sum(T[Tcenter-dr:Tcenter+dr,Tcenter-dr:Tcenter+dr])
            for mt in tqu_template_maps[band].maps:
                tqu_template_maps[band].maps[mt] *= (fluxt/fluxT)
            print(band,fluxt/fluxT)
    if debug:
        bandfig = {'90':1,'150':2,'220':3}
        plt.figure(figsize=(15,20))
        for band in perband_T_coadds:
            t = perband_T_coadds[band]
            
            T = tqu_template_maps[band].maps['T']
            
            ## Uncomment to interpolate and align maps
            xlen = 2*tcenter
            datax,datay  = np.meshgrid(np.arange(xlen),np.arange(xlen))
            gridpoints  = np.linspace(0,xlen,240)
            gridx,gridy = np.meshgrid(gridpoints,gridpoints)
            t = griddata((datax.flatten(),datay.flatten()),perband_T_coadds[band].flatten(),(gridx,gridy))
            T = griddata((datax.flatten(),datay.flatten()),tqu_template_maps[band].maps['T'].flatten(),(gridx,gridy))
            T = optimize_map_alignment(T,t,band,5,5)
            
            plt.subplot(4,3,bandfig[band])
            plt.imshow(T)
            plt.colorbar(fraction=.046,pad=.04)
            plt.title('%sGHz normed TQU templates'%(band))
            plt.grid(color='w',ls='--')

            plt.subplot(4,3,3+bandfig[band])
            plt.imshow(t)
            plt.colorbar(fraction=.046,pad=.04)
            plt.title('%sGHz T Coadds '%(band))
            plt.grid(color='w',ls='--')

            plt.subplot(4,3,6+bandfig[band])
            plt.imshow(T-t)
            plt.colorbar(fraction=.046,pad=.04)
            plt.title('%sGHz  normed TQU template - T Coadd '%(band))
            plt.grid(color='w',ls='--')
            
            avg = (T+t)/2.
            frac_res = (T-t)/avg
            frac_res[avg<1.] = 0.
            plt.subplot(4,3,9+bandfig[band])
            plt.imshow(100.*frac_res,vmin=-10,vmax=10)
            plt.colorbar(fraction=.046,pad=.04)
            plt.title('%sGHz  residual percent of avg '%(band))
            plt.grid(color='w',ls='--')
        plt.suptitle('%s Per-Band Coadd Residuals %s'%(source,obs))
        plt.savefig(savedir+'%s_PerBand_Coaddmap_Residuals_%s.png'%(source,obs))
        plt.close('all')
        
    return tqu_template_maps
    
    
    

def get_parameter_range(minval,valrange,n = 11):
    '''
    Used to change the parameter range for minimizing chi squared fits.
    Takes in the previous min and the previous range; then outputs a
    linear interpolation between the nearest neighbors of the previous range.
    
    Arguments
    ---------
    minval : float
        Parameter value which previously minimized chi squared.
    
    valrange : np array
        Parameter values used in previous chi squared fit.
        
    n : int
        Number of elements to use in the interpolation.  
     
        
    Returns
    -------
    output name : output type
        explanation of output
    '''
    if len(valrange)>1:
            dval = np.mean(np.diff(valrange))
            valrange = np.linspace(minval-dval,minval+dval,n)
    return valrange


def calc_chi2(grange,prange,thetarange,t,T,Q,U,var,dof):
    '''
    Returns a 3D numpy array of the chi squared landscape for the input parameters.
    Also returns the minimized chi2 of this landscape.
    
    
    Arguments
    ---------
    Grange : list or array
        List of gain values which are being used in the fit.
        
    prange : list or array
        List of polarization efficiency values which are being used in the fit.
        
    Thetarange : list or array
        List of polarization angle values which are being used in the fit.
        
    t : 2D numpy array
        Temperature map to be fit to template TQU
     
    T,Q,U : 2D numpy arrays
        Template Maps (typically per-band coadds) to use in fitting.
        
    var : float
        Map variance, for use in the calculation of chi2
       
    dof : int
        Number of degrees of freedom (number of pixels in the fit)
        
        
    Returns
    -------
    chi2 : 3D numpy array
        Chi squared landscape of the input parameter matrix.
    
    '''
    chi2 = np.zeros((len(grange),len(prange),len(thetarange)))
    i = 0
    for G in grange:
        j = 0
        for p in prange:
            k = 0
            for theta in thetarange:
                model_t = np.asarray(.5* G * ((2-p)*T + p * np.cos(2*theta) * Q + p * np.sin(2*theta) * U))
                chi2[i][j][k] = np.nansum((t - model_t)**2/var)
                k+=1
            j+=1
        i+=1
    
    return chi2




def likelihood_prior(x,mu,sigma):
    mu = np.deg2rad(mu)
    sigma=np.deg2rad(sigma)
    N = (1./(sigma*np.sqrt(2*np.pi))) * np.exp((-(CF.difference(x,mu))**2)/(2*sigma**2))
    return N
    
def fit_polarization_params(nomang,band,t,T,Q,U,mask,thetamean,dtheta=1.,gmean=1.5,dg=.1,pmean=1.0,dp=.1,VARIABLE = 'gain'):
    '''
    Finds the parameters which minimize chi2 fit of t to template TQU maps within 
    the region masked by mask.
    
    One can choose to marginalize Theta over detector gain, or polarization efficiency,
    or both.
    
    Can plot the chi2 landscape, log likelihood landscape and marginalized distributions
    by choosing to set plotname. One should choose a full path for the plot name, 
       for example:  plotname = '/path/to/plot/directory/90GHz_1degNomang' 
    and this function will append the plot name to the end of this path.
    
    Arguments
    ---------
    band : str
        Band of TQU maps; i.e. '90', '150', or '220'
        
    t : 2D numpy array
        Temperature map to be fit to template TQU
     
    T,Q,U : 2D numpy arrays
        Template Maps (typically per-band coadds) to use in fitting.    
        
    mask : 2D numpy array
        Mask array for masking the region you want to fit.
        
    VARIABLE : str
        Variable which is to be fit; either 'G','p' or 'both'
        for gain, pol eff, and both. This variable gets marginalized
        over during the calculation of best fit theta and uncerts.
        
    plotname : str
        Path to and plot name for plotting the chi2 landscape, chi2 hist
        and log-likelihood landscape.
        
    Returns
    -------
    Marginalized parameters : list
        List of the mean marginalized best-fit parameters and their standard deviations.
        returned as [[gain,gain_std],[p,p_std],[theta,theta_std]]
    '''
    # Bolo Map t given by: 
    # t = 1/2 G * [(2-p)*T + p * cos(2Theta) * Q + p * sin(2Theta) * U]
    # G is gain, p is pol eff, Theta is pol angle of bolometer
    
    ## Generate inverse mask for generating noise variance.
    inverse_mask = np.copy(1.-mask)
    dx = np.shape(inverse_mask)[0]//2

    #make sure nucleus is not included in noise calc (if nucleus masked out)
    inverse_mask[dx-9:dx+9,dx-9:dx+9] *= 0
    flatmask = np.ndarray.flatten(mask)
    inverse_mask = np.ndarray.flatten(inverse_mask)
    dof = len(np.nonzero(flatmask)[0])
    #var = np.nanvar(t[50:,:10])+np.nanvar(T[50:,:10])+np.nanvar(Q[50:,:10])+np.nanvar(U[50:,:10])
    #print('corner var: %.3f'%var)
    ## apply source mask to template maps and t map.
    T = np.ndarray.flatten(T)*flatmask
    Q = np.ndarray.flatten(Q)*flatmask
    U = np.ndarray.flatten(U)*flatmask
    t = np.ndarray.flatten(t)
    nonzeronoise = t[np.nonzero(inverse_mask)[0]]
    ## Get sum variance of maps being fit. used for fit; assuming gaussian random noise
    var = np.nanvar(nonzeronoise)+np.nanvar(T[np.nonzero(inverse_mask)[0]])\
        +np.nanvar(Q[np.nonzero(inverse_mask)[0]]) + np.nanvar(U[np.nonzero(inverse_mask)[0]])
    
    
    t *= flatmask
    
    
    ## Ideal cases; used to set one while fitting the other.
    gband = {'90':2.,'150':2.,'220':2.}
    pband = {'90':1.,'150':1.,'220':1.}
    dthetaband = {'90':30.,'150':30.,'220':30.}
    prior_width = {'90':60.,'150':60.,'220':60.}
    nbins = 6
   
    
    if VARIABLE == 'both':
        variables = ['gain','rho']
    elif VARIABLE == 'fixed':
        variables = ['fixed']
    else:
        variables = [VARIABLE]
    
    dof-=(1+len(variables))
    
    
    if VARIABLE == 'both':
        thetarange = np.deg2rad(np.arange(thetamean-nbins*dtheta,thetamean+nbins*dtheta,dtheta))
        grange = np.arange(max(gmean-nbins*dg,0),gmean+nbins*dg,dg)
        prange = np.arange(max(pmean-nbins*dp,0),min(pmean+nbins*dp,1),dp)
        #print(nomang)
        #print(np.rad2deg(thetarange))
        #print(grange)
        #print(prange)
        
        chi2 = calc_chi2(grange,prange,thetarange,t,T,Q,U,var,dof)
        #print('min red chi2: %s'%(np.amin(chi2)/dof))
        minchi2=np.where(chi2==np.amin(chi2))
        #print(minchi2[0])
        return [[grange[minchi2[0][0]],0],[prange[minchi2[1][0]],0],[np.rad2deg(thetarange[minchi2[2][0]]),0]]
    
    ## Marginalize Theta over the variables.
    ## Gain, pol eff have range of +- .5
    for variable in variables:
        if variable == 'fixed':
            grange = [gband[band]]
            prange = [pband[band]]
            
        elif variable == 'gain' or variable == 'G' or variable == 'g':
            grange = np.linspace(gmean-nbins*dg,gmean+nbins*dg,dg)
            prange = [pband[band]]
        else:
            prange = np.linspace(pmean-nbins*dp,pmean+nbins*dp,dp)
            grange = [gband[band]]

        thetarange = np.deg2rad(np.arange(thetamean-nbins*dtheta,thetamean+nbins*dtheta,dtheta))
        

        ## Calc Chi^2 
        
        chi2 = calc_chi2(grange,prange,thetarange,t,T,Q,U,var,dof)
        #print('min red chi2: %s'%(np.amin(chi2)/dof))
        minchi2=np.where(chi2==np.amin(chi2))
        #print(minchi2[0])
        return [[grange[minchi2[0][0]],0],[prange[minchi2[1][0]],0],[np.rad2deg(thetarange[minchi2[2][0]]),0]]
    
       
    
    
    
def fit_polarization_params_old(nomang,band,t,T,Q,U,mask,VARIABLE = 'gain',plotname=''):
    '''
    Finds the parameters which minimize chi2 fit of t to template TQU maps within 
    the region masked by mask.
    
    One can choose to marginalize Theta over detector gain, or polarization efficiency,
    or both.
    
    Can plot the chi2 landscape, log likelihood landscape and marginalized distributions
    by choosing to set plotname. One should choose a full path for the plot name, 
       for example:  plotname = '/path/to/plot/directory/90GHz_1degNomang' 
    and this function will append the plot name to the end of this path.
    
    Arguments
    ---------
    band : str
        Band of TQU maps; i.e. '90', '150', or '220'
        
    t : 2D numpy array
        Temperature map to be fit to template TQU
     
    T,Q,U : 2D numpy arrays
        Template Maps (typically per-band coadds) to use in fitting.    
        
    mask : 2D numpy array
        Mask array for masking the region you want to fit.
        
    VARIABLE : str
        Variable which is to be fit; either 'G','p' or 'both'
        for gain, pol eff, and both. This variable gets marginalized
        over during the calculation of best fit theta and uncerts.
        
    plotname : str
        Path to and plot name for plotting the chi2 landscape, chi2 hist
        and log-likelihood landscape.
        
    Returns
    -------
    Marginalized parameters : list
        List of the mean marginalized best-fit parameters and their standard deviations.
        returned as [[gain,gain_std],[p,p_std],[theta,theta_std]]
    '''
    # Bolo Map t given by: 
    # t = 1/2 G * [(2-p)*T + p * cos(2Theta) * Q + p * sin(2Theta) * U]
    # G is gain, p is pol eff, Theta is pol angle of bolometer
    
    ## Generate inverse mask for generating noise variance.
    inverse_mask = np.copy(1.-mask)
    dx = np.shape(inverse_mask)[0]//2

    #make sure nucleus is not included in noise calc (if nucleus masked out)
    inverse_mask[dx-9:dx+9,dx-9:dx+9] *= 0
    flatmask = np.ndarray.flatten(mask)
    inverse_mask = np.ndarray.flatten(inverse_mask)
    dof = len(np.nonzero(flatmask)[0])
    #var = np.nanvar(t[50:,:10])+np.nanvar(T[50:,:10])+np.nanvar(Q[50:,:10])+np.nanvar(U[50:,:10])
    #print('corner var: %.3f'%var)
    ## apply source mask to template maps and t map.
    T = np.ndarray.flatten(T)*flatmask
    Q = np.ndarray.flatten(Q)*flatmask
    U = np.ndarray.flatten(U)*flatmask
    t = np.ndarray.flatten(t)
    nonzeronoise = t[np.nonzero(inverse_mask)[0]]
    ## Get sum variance of maps being fit. used for fit; assuming gaussian random noise
    var = np.nanvar(nonzeronoise)+np.nanvar(T[np.nonzero(inverse_mask)[0]])\
        +np.nanvar(Q[np.nonzero(inverse_mask)[0]]) + np.nanvar(U[np.nonzero(inverse_mask)[0]])
    
    
    t *= flatmask
    
    
    ## Ideal cases; used to set one while fitting the other.
    gband = {'90':2.,'150':2.,'220':2.}
    pband = {'90':1.,'150':1.,'220':1.}
    dthetaband = {'90':30.,'150':30.,'220':30.}
    prior_width = {'90':60.,'150':60.,'220':60.}
    ## Set nominal output values
    gmean = gband[band]
    gstd  = 0.
    pmean = pband[band]
    pstd  = 0.
    
    angstart = nomang-dthetaband[band]  
    angend   = nomang+dthetaband[band]

    ## resolution of 2 deg
    ## UPDATE: MANUALLY ENTERED NBINS BELOW!!!
    #nbins    = 61

    marginalized_theta = {}
    
    if VARIABLE == 'both':
        variables = ['gain','rho']
    elif VARIABLE == 'fixed':
        variables = ['fixed']
    else:
        variables = [VARIABLE]
    
    dof-=(1+len(variables))
    
    
    if VARIABLE == 'both':
        thetarange = np.deg2rad(np.linspace(angstart,angend,121))
        grange = np.linspace(gband[band]-.3,gband[band]+.3,61)
        prange = np.linspace(pband[band] - .3,pband[band]+.3,61)    
        chi2 = calc_chi2(grange,prange,thetarange,t,T,Q,U,var,dof)
        #print('min red chi2: %s'%(np.amin(chi2)/dof))
        minchi2=np.where(chi2==np.amin(chi2))
        #print(minchi2[0])
        return [[grange[minchi2[0][0]],0],[prange[minchi2[1][0]],0],[np.rad2deg(thetarange[minchi2[2][0]]),0]]
    
    ## Marginalize Theta over the variables.
    ## Gain, pol eff have range of +- .5
    for variable in variables:
        if variable == 'fixed':
            grange = [gband[band]]
            prange = [pband[band]]
            
        elif variable == 'gain' or variable == 'G' or variable == 'g':
            grange = np.linspace(gband[band]-.3,gband[band]+.3,61)
            prange = [pband[band]]
        else:
            prange = np.linspace(pband[band] - .3,pband[band]+.3,61)
            grange = [gband[band]]

        thetarange = np.deg2rad(np.linspace(angstart,angend,181))


        ## Calc Chi^2 
        
        chi2 = calc_chi2(grange,prange,thetarange,t,T,Q,U,var,dof)
        #print('min red chi2: %s'%(np.amin(chi2)/dof))
        minchi2=np.where(chi2==np.amin(chi2))
        #print(minchi2[0])
        return [[grange[minchi2[0][0]],0],[prange[minchi2[1][0]],0],[np.rad2deg(thetarange[minchi2[2][0]]),0]]
    
        #print(nomang)
        
        ## set Theta as x and param to marginlize over as y
        ## Set 2D residuals space --> res[G][p][theta]
        x = thetarange

        if variable == 'gain' or variable == 'G' or variable == 'g':

            y = grange
            res2d = chi2[:,0,:]
            ylabel = 'Gain'
        else:
            y = prange
            res2d = chi2[0,:,:]
            ylabel = 'Pol eff'
            
        
        ## Get Likelihoods and marginalized 1D arrays
        '''
        print('dof: %s'%dof)
        print('mean,var from scipy: %.3f, %.3f'%(scipy_chi2.stats(dof,moments='mv')))
        
        residualfit = t-.5*grange[minchi2[0][0]]*prange[minchi2[1][0]]*(T*(2-prange[minchi2[1][0]])/prange[minchi2[1][0]] + np.cos(2*theta) * Q + np.sin(2*theta) * U)
        print('mapvar: %s'%var)
        scale   = np.amin(res2d)
        
        minchi2=np.where(chi2==np.amin(chi2))
        print(band,nomang)
        theta = np.rad2deg(thetarange[minchi2[2][0]])
        print('min chi2 theta: %.2f'%theta)
        '''
        #print('%sGHz %s deg, %.2f'%(band,nomang,np.amin(res2d)/dof))
        likelihood = scipy_chi2.pdf(res2d,dof,scale=np.amin(res2d/dof))# np.exp(-res2d/2) # 
        likelihood/=np.nansum(likelihood)
        
        
        marginalized_x = np.sum(likelihood,axis = 0) #* likelihood_prior(x,nomang,prior_width[band])
        marginalized_y = np.sum(likelihood,axis = 1)
        ## save for marginalizing over both p and gain
        marginalized_theta[variable]= np.copy(marginalized_x)
        
        ## Finessing angles to get weighted mean,variance 
        #x*=2
        margmeanx = CF.nanmean(x,w=marginalized_x)#/2
        margstdx = np.sqrt(CF.nanvar(x,w=marginalized_x))#/2
        #x/=2
        
        margmeany = np.sum(marginalized_y*y)
        margstdy  = np.sqrt(np.sum(marginalized_y*(y-margmeany)**2))
        

        if variable == 'G' or variable == 'gain' or variable =='g':
            gmean = margmeany
            gstd  = margstdy
        else:
            pmean = margmeany
            pstd  = margstdy
    ## Marginalizing over both variables
    if VARIABLE =='both':
        comarg_x = np.ones(len(marginalized_x))
        for variable in marginalized_theta:
            comarg_x *= marginalized_theta[variable]**2
        comarg_x = np.sqrt(comarg_x)
        
        margmeanx = CF.nanmean(x,w=comarg_x)
        margstdx = np.sqrt(CF.nanvar(x,w=comarg_x))
        
    else:
        comarg_x = marginalized_x
    comarg_x/=np.sum(comarg_x)
    '''
    if margmeanx<0.:
        margmeanx+=np.pi
    '''
    ## Useful for debugging; if you want to plot, set plotname to 
    ## something useful that explains what bolos you're plotting.
    if plotname!='':
        if variable == 'fixed':
            plt.figure(figsize=(10,5))
            #y = np.log10(likelihood)[0]
            plt.plot(likelihood[0])
            ## Get upper and lower 1-sigma from mean.
            lowersigma=0
            uppersigma=0
            integral = 0
            mux = 0
            integral = 0
            i = 0
            while integral < .5:
                mux = i
                integral+= likelihood[0][i]
                i+=1
            i =0
            integral = 0
            dl = np.diff(likelihood[0])[0]
            while integral < .16:
                integral+=likelihood[0][i]
                lowersigma = i
                i+=1
            integral = 0
            while integral < .68:
                integral+=likelihood[0][i]
                uppersigma = i
                i+=1
            plt.plot([lowersigma,lowersigma],[0,max(likelihood[0])],color='r')
            plt.plot([uppersigma,uppersigma],[0,max(likelihood[0])],color='r')
            mu = x[mux]
            psigma = np.rad2deg(abs(CF.difference(x[uppersigma],mu)))
            msigma = np.rad2deg(abs(CF.difference(x[lowersigma],mu)))
            mu = np.rad2deg(mu)
            plt.plot([mux,mux],[0,max(likelihood[0])],color='b',label=r'%.2f$^{+%.2f}_{-%.2f}$'%(mu,psigma,msigma))
            #print(y)
            #plt.plot(res2d[0],label='chi2')
            #plt.yticks(range(0,len(y),1),np.around(y,decimals=1),rotation=0,fontsize=6)
            plt.xticks(np.arange(0,len(x)),np.around(np.rad2deg(x),decimals=1),rotation=90,fontsize=6)
            plt.grid()
            plt.ylabel('Probability')
            plt.xlabel('Angle (deg)')
            plt.legend()
            plt.title(plotname)
            plt.savefig('%s_loglikelihood.png'%plotname)
            plt.close('all')
            ## Plot Chi^2 vs angle
            plt.figure(figsize=(10,5))
            y = chi2[0,0,:]/dof
            plt.plot(y)
            minchi2 = np.min(y)
            sigmaline = minchi2+2.3/dof
            plt.plot([0,len(y)],[sigmaline,sigmaline])
            
            #plt.plot(res2d[0],label='chi2')
            #plt.yticks(range(0,len(y),1),np.around(y,decimals=1),rotation=0,fontsize=6)
            plt.xticks(np.arange(0,len(x)),np.around(np.rad2deg(x),decimals=1),rotation=90,fontsize=6)
            plt.grid()
            plt.ylabel(r'$\chi^2 / \nu$')
            plt.xlabel('Angle (deg)')
            plt.title(plotname)
            plt.savefig('%s_redChi2.png'%plotname)
            plt.close('all')
            
        else:
            ## Plot histogram of flattened chi^2
            plt.figure()
            plt.hist(np.ndarray.flatten(res2d/dof),bins=251,histtype='step',color='k')
            plt.grid()
            plt.yscale('log')
            plt.xscale('log')
            plt.xlabel('Reduced $\chi^2$')
            plt.title(plotname)
            plt.savefig('%s_Chi2_Hist.png'%plotname)
            plt.close('all')


            ## If num params = 2 ; deltaX^2 = [2.3,4.61,9.21]
            ##               = 3 ; [3.5,6.25,11.3]
            ## Plot chi^2 landscape with 1,2,3 sigma contours
            plt.figure(figsize=(20,10))
            plt.imshow(res2d/dof,cmap='bone',vmin=np.amin(res2d/dof),vmax=np.amin(res2d/dof)+5*np.sqrt(2/dof))
            plt.colorbar(fraction=.046,pad=.04,label='Reduced $\chi^2$')
            plt.yticks(range(0,len(x),1),np.around(y,decimals=3),rotation=0,fontsize=6)
            plt.xticks(range(0,len(x),1),np.around(np.rad2deg(x),decimals=1),rotation=90,fontsize=6)
            plt.ylabel(ylabel)
            plt.xlabel('Angle (deg)')
            if VARIABLE == 'p' or VARIABLE=='both':
                plt.scatter(minchi2[2][0],minchi2[1][0],marker='x',color='w',label='MinChi Theta: %.1f deg'%np.rad2deg(thetarange[minchi2[2][0]]))
                
            else:
                plt.scatter(minchi2[2][0],minchi2[0][0],marker='x',color='w',label='MinChi Theta: %.1f deg'%np.rad2deg(thetarange[minchi2[2][0]]))
            levels = np.amin(res2d/dof)+np.array([1.,2.,3.])/(np.sqrt(2*dof)) #+np.array([2.3,4.61,9.21])*(1./dof)
            cset = plt.contour(res2d/dof,levels,colors='w')
            clabelnames = {}
            for i in range(len(levels)):
                clabelnames[levels[i]] = '%i-$\sigma$'%(i+1)
            plt.clabel(cset,fmt=clabelnames,fontsize=12)
            plt.title('Reduced $\chi^2$ Landscape %s'%plotname.split('.')[0])
            plt.legend()
            plt.savefig('%s_Chi2_Landscape.png'%plotname)
            plt.close('all')

            ## Plot log-likelihood landscape.
            plt.figure(figsize=(20,10))
            plt.imshow(likelihood,cmap='bone_r')
            plt.colorbar(fraction=.046,pad=.04,label='Probability')
            xlen = len(x)
            ylen = len(y)
            ## Get upper and lower 1-sigma from mean.
            lowersigma=0
            uppersigma=0
            integral = 0
            mux = 0
            integral = 0
            i = 0
            if VARIABLE=='both':
                margval=comarg_x
            else:
                margval=marginalized_x
            while integral < .5:
                mux = i
                integral+= margval[i]
                i+=1
            i =0
            integral = 0
            while integral < .16:
                integral+=margval[i]
                lowersigma = i
                i+=1
            integral = 0
            while integral < .68:
                integral+=margval[i]
                uppersigma = i
                i+=1
            plt.plot([lowersigma,lowersigma],[0,ylen],color='r')
            plt.plot([uppersigma,uppersigma],[0,ylen],color='r')
            mu = x[mux]
            psigma = np.rad2deg(abs(CF.difference(x[uppersigma],mu)))
            msigma = np.rad2deg(abs(CF.difference(x[lowersigma],mu)))
            mu = np.rad2deg(mu)
            plt.plot([mux,mux],[0,ylen],color='g',label=r'%.2f$^{+%.2f}_{-%.2f}$'%(mu,psigma,msigma))
            
            ## Will plot co-marginalized lines over top of the log-likelihood plot
            if VARIABLE == 'both':
                plt.plot(range(0,len(x),1),(1.-comarg_x/np.max(comarg_x))*ylen,color='r',label='co-Marginalized Theta: %4.2f +- %4.2f'%(np.rad2deg(margmeanx),np.rad2deg(margstdx)) )
            else:
                plt.plot(range(0,len(x),1),(1.-marginalized_x/np.max(marginalized_x))*ylen,color='r',label='Marginalized Theta: %4.2f +- %4.2f'%(np.rad2deg(margmeanx),np.rad2deg(margstdx)) )
            plt.plot(marginalized_y/np.max(marginalized_y)*xlen/2,range(0,len(y),1),color='b',label='Marginalized %s: %4.2f +- %4.2f'%(ylabel,margmeany,margstdy) )
            plt.legend()
            plt.yticks(range(0,len(y),1),np.around(y,decimals=3),rotation=0,fontsize=6)
            plt.xticks(range(0,len(x),1),np.around(np.rad2deg(x),decimals=1),rotation=90,fontsize=6)
            plt.ylabel(ylabel)
            plt.xlabel('Angle (deg)')
            #plt.scatter(minvals[2][0],minvals[0][0],marker='x',color='w')
            '''
            levels = np.amin(res2d)+np.array([2.3,4.61,9.21])
            cset = plt.contour(res2d,levels,colors='w')
            clabelnames = {}
            for i in range(len(levels)):
                clabelnames[levels[i]] = '%i-$\sigma$'%(i+1)
            plt.clabel(cset,fmt=clabelnames,fontsize=12)
            '''
            plt.title('Probability Landscape %s'%plotname.split('.')[0])
            plt.savefig('%s_Likelihood_Landscape.png'%plotname)
            plt.close('all')
        
    ## return Gain,pol efficiency,pol angle  
    #print([gmean,gstd],[pmean,pstd],[np.rad2deg(margmeanx),np.rad2deg(margstdx)])
    return [[gmean,gstd],[pmean,pstd],[np.rad2deg(margmeanx),np.rad2deg(margstdx)]]
            

def fit_individual_bolometer_polarization(bolo_maps,template_tqu_maps,mask,log,goodbolos = [],mapsize = 30,variable='gain',plot_loc='./tmp_plots/',PLOT=0,rot_ang = 0.):
    '''
    This function fits single-bolo maps to the template TQU maps.
    
    Arguments
    ---------
    arg name : arg type
        explanation of arg
        
    Returns
    -------
    output name : output type
        explanation of output
    '''
    
    ## Fits nominal angle maps per wafer, per band to TQU sum maps.  
    BANDS = ['90','150','220']
    
    t = {}
    q = {}
    u = {}
    for band in BANDS:
        if band in template_tqu_maps:
            coband = band
        elif int(band) in template_tqu_maps:
            coband = int(band)
        else:
            coband = band +'GHz'
        summapsize = np.shape(template_tqu_maps[coband].maps['T'])[0]//2
        t[band] = template_tqu_maps[coband].maps['T'] # perbandcoadds[band] #
        q[band] = template_tqu_maps[coband].maps['Q']
        u[band] = template_tqu_maps[coband].maps['U']
        
    if len(goodbolos) == 0:
        goodbolos = bolo_maps.keys()
        log.write('\n !!!!!!!!!!!!! \n')
        log.write('No Good Bolos Provided, using all bolos in bolo_maps\n !!!!!!!!!!!! \n')
        
    for band in mask:
        ms = np.shape(mask[band])[0]//2
        ## if mask is smaller than set map size; will reduce map size.
        if ms <= mapsize:
            mapsize=ms
        else:
            mask[band] = mask[band][ms-mapsize:ms+mapsize,ms-mapsize:ms+mapsize]
    
    numbolos = len(goodbolos)
    startoffit = time.time()
    bolonum = 0.
    percentdone = 1.
    plotnum = 0
    fitdata = {}
    for bolo in goodbolos:  
        bolonum+=1.
        if (bolonum/numbolos)*100.>percentdone:
            dt = (time.time()-startoffit)
            log.write('\n ------------------------------------------------------\n')
            log.write('%i Percent Done!\n'%percentdone)
            log.write('Elapsed time in calc: %i s.\n'%dt)
            log.write('Average time per bolo: %.3f s.\n'%((dt/float(bolonum))))
            log.write('Estimated time remaining: %.3f s.\n'%(dt/float(bolonum)*(numbolos-bolonum)))
            log.write(' ------------------------------------------------------\n')
            percentdone+=1.
                  
        fitdata[bolo] = {}
        ang  = bolo_maps[bolo].nomang
        band = bolo_maps[bolo].band
        T = t[band][summapsize-mapsize:summapsize+mapsize,summapsize-mapsize:summapsize+mapsize]
        Q = q[band][summapsize-mapsize:summapsize+mapsize,summapsize-mapsize:summapsize+mapsize]
        U = u[band][summapsize-mapsize:summapsize+mapsize,summapsize-mapsize:summapsize+mapsize]
        
        cmapsize = bolo_maps[bolo].shape
        if cmapsize[0]//2>= mapsize and cmapsize[1]//2 >= mapsize :
            xsize = cmapsize[0]//2
            ysize = cmapsize[1]//2
            ti = bolo_maps[bolo].maps['T'][xsize-mapsize:xsize+mapsize,ysize-mapsize:ysize+mapsize]
            ## Scale individual maps to the mapsize
            if np.all(~np.isfinite(ti)):
                        continue
            ti = np.nan_to_num(ti)
        
            ## Do the fitting
            if PLOT:
                plotname = plot_loc+'%s.png'%bolo
            else:
                plotname = ''
                
            nomang = ang+rot_ang
            
            
            dtheta=5.
            dg = .05
            dp = .05
            p2=fit_polarization_params(nomang,band,ti,T,Q,U,mask[band],thetamean=nomang,dtheta=dtheta,
                                       gmean=2.0,dg=dg,pmean=.8,dp=dp,VARIABLE = variable)
            i=0
            #print('0th fit: ',p2)
            while True:
                res = CF.difference(p2[2][0],nomang,deg=True)
                res0 = np.copy(res)
                newnomang = nomang
                while abs(res)>=30.:
                    i+=1
                    last_res = np.copy(res)
                    p2=fit_polarization_params(nomang,band,ti,T,Q,U,mask[band],thetamean=p2[2][0],dtheta=dtheta,
                                       gmean=p2[0][0],dg=dg,pmean=p2[1][0],dp=dp,VARIABLE = variable)
                    res = CF.difference(p2[2][0],nomang,deg=True)
                    #print('%s iteration: '%i,p2)
                    if abs(last_res)<=abs(res):
                        res = last_res
                        break
                    if i > 10:
                        break
                if dtheta>1.:
                    i+=1
                    dtheta/=2.
                    dg /= 2.
                    dp /= 2.
                    p2=fit_polarization_params(nomang,band,ti,T,Q,U,mask[band],thetamean=p2[2][0],dtheta=dtheta,
                                       gmean=p2[0][0],dg=dg,pmean=p2[1][0],dp=dp,VARIABLE = variable)
                    res = CF.difference(p2[2][0],nomang,deg=True)
                    #print('%s iteration: '%i,p2)
                else:
                    break
            #print('best fit: ',p2)         
            
            
            ## Rotate the angle into [0,180.]
            ## Rotate phi by 90 degrees since 5/30/2019 HWM update
            ## included a fit to optimize rotation angle
            
            
            p2[2][0] -= rot_ang
            
            if p2[2][0]<0:
                p2[2][0]+=180.
            if p2[2][0]>180.:
                p2[2][0]-=180.
            
            
            
            ## For debugging; manually change this
            if plotnum < 0:
                    plt.figure(figsize=(40,8))
                    plt.subplot(151)
                    plt.imshow(ti*mask[str(band)])
                    plt.title('%s t_i'%bolo)
                    plt.colorbar()
                    plt.grid(color='w',ls='--')
                    plt.subplot(152)
                    plt.imshow(T*mask[str(band)])
                    plt.colorbar()
                    plt.grid(color='w',ls='--')
                    plt.title('T')
                    plt.subplot(153)
                    plt.imshow((ti-T)*mask[str(band)])
                    plt.colorbar()
                    plt.grid(color='w',ls='--')
                    plt.title('t_i-T')
                    plt.subplot(154)
                    plt.imshow((.5* p2[0] * ((2-p2[1])*T + p2[1] * np.cos(2*np.deg2rad(p2[2])) * Q + p2[1] * np.sin(2*np.deg2rad(p2[2])) * U))*mask[str(band)])
                    plt.title("Best Fit t_i")
                    plt.colorbar()
                    plt.grid(color='w',ls='--')
                    plt.subplot(155)
                    residuals = (ti-.5* p2[0] * ((2-p2[1])*T + p2[1] * np.cos(2*np.deg2rad(p2[2])) * Q +
p2[1] * np.sin(2*np.deg2rad(p2[2])) * U))
                    var = np.nanvar(ti[np.nonzero(ti*mask[str(band)])])
                    plt.hist(np.ndarray.flatten(residuals)/var,bins = 51)
                    plt.title("Fit Residual / var")
                    plt.grid(color='w',ls='--')
                    plt.savefig(plot_loc+'%s_FittingPlot.png'%bolo)
                    plt.close('all')
                    plotnum+=1
                    
            phi    = p2[2][0]
            nomang =  ang 
            residual_angle = CF.difference(phi,nomang,deg=True)

            fitdata[bolo]['G']      = p2[0]
            fitdata[bolo]['rho']    = p2[1]
            fitdata[bolo]['rawPhi'] = p2[2]
            fitdata[bolo]['nom']    = nomang
            #fitdata[bolo]['newnom'] = newnomang
            fitdata[bolo]['rotPhi'] = phi
            fitdata[bolo]['res']    = residual_angle
            
        

    return fitdata





def fit_nominal_angle_coadd_polarization(nomangle_maps,perbandcoadds,template_tqu_maps,mask,numbolos,log,mapsize=30,variable = 'p',plot_loc='./tmp_plots/',PLOT = 0,rot_ang = 0.):

    '''
    This function fits maps made from per-wafer, per-band, per-nominal angle coadds to the template TQU maps.
    
    Arguments
    ---------
    arg name : arg type
        explanation of arg
        
    Returns
    -------
    output name : output type
        explanation of output
    '''
    ## Fits nominal angle maps per wafer, per band to band-averaged TQU template maps.  
    BANDS = ['90','150','220']
    fitdata = {}
    t = {}
    q = {}
    u = {}
    for band in mask:
        ms = np.shape(mask[band])[0]//2
        if ms <= mapsize:
            mapsize=ms
        else:
            mask[band] = mask[band][ms-mapsize:ms+mapsize,ms-mapsize:ms+mapsize]
    for band in BANDS:
        if band in template_tqu_maps:
            coband = band
        elif int(band) in template_tqu_maps:
            coband = int(band)
        else:
            coband = band +'GHz'
        cmapsize = np.shape(template_tqu_maps[coband].maps['T'])[0]//2
        t[band] = template_tqu_maps[coband].maps['T'] # perbandcoadds[band] #
        q[band] = template_tqu_maps[coband].maps['Q']
        u[band] = template_tqu_maps[coband].maps['U']    
    for wafer in sorted(nomangle_maps.keys())[::-1]:    
        log.write('Fitting Nominal Angle Coadds for %s\n'%wafer)
        fitdata[wafer] = {}
        for ang in sorted(nomangle_maps[wafer].keys()):
            fitdata[wafer][ang] = {}
            for band in nomangle_maps[wafer][ang]:
                if numbolos[wafer][ang][band]>0:
                    band=str(band)
                    fitdata[wafer][ang][band] = {}
                    cmapsize = np.shape(nomangle_maps[wafer][ang][band])[0]//2
                    log.write('%s %sGHz %i deg \n'%(wafer,band,ang))
                    log.write('Number of bolos: %i\n'%numbolos[wafer][ang][band])
                    log.write('Nominal Phi = %.2f\n'%(ang))
                    
                    ## Make mask and maps same size
                    if cmapsize<mapsize or cmapsize<ms:
                        mapsize = cmapsize
                        shrunkmask = mask[band][ms-cmapsize:ms+cmapsize,ms-cmapsize:ms+cmapsize]
                    else:
                        shrunkmask = mask[band]
                    ## Make sure all maps are equivalently shaped.
                    T = t[band][cmapsize-mapsize:cmapsize+mapsize,cmapsize-mapsize:cmapsize+mapsize]
                    Q = q[band][cmapsize-mapsize:cmapsize+mapsize,cmapsize-mapsize:cmapsize+mapsize]
                    U = u[band][cmapsize-mapsize:cmapsize+mapsize,cmapsize-mapsize:cmapsize+mapsize]
                    
                    ti = nomangle_maps[wafer][ang][band][cmapsize-mapsize:cmapsize+mapsize,cmapsize-mapsize:cmapsize+mapsize]
                    ## If there are NaNs in the template maps; something went wrong. 
                    if np.all(~np.isfinite(ti)) or np.all(~np.isfinite(T)) or np.all(~np.isfinite(Q)) or np.all(~np.isfinite(U)):
                        continue
                        
                    
                    
                    
                    nomang = ang+rot_ang
                    
                    
                    if PLOT:
                        plotname = plot_loc+'%s_%sGHz_%sdeg_%s'%(wafer,band,nomang,variable)
                    else:
                        plotname = ''
                    
                    
                    
                    
                    dtheta=3.
                    dg = .05
                    dp = .05
                    
                    
                    p2=fit_polarization_params(nomang,band,ti,T,Q,U,mask[band],thetamean=nomang,dtheta=dtheta,
                                       gmean=2.0,dg=dg,pmean=.9,dp=dp,VARIABLE = variable)
                    
                    '''
                    plt.figure(figsize=(40,16))
                    plt.subplot(241)
                    plt.imshow(ti*mask[band])
                    plt.title('%s t_i'%plotname)
                    plt.colorbar(fraction=.046,pad=.04)
                    plt.grid(color='w',ls='--')
                    plt.subplot(242)
                    plt.imshow(T*mask[band])
                    plt.colorbar(fraction=.046,pad=.04)
                    plt.grid(color='w',ls='--')
                    plt.title('T')
                    plt.subplot(243)
                    plt.imshow((ti-.5* p2[0][0] * ((2-p2[1][0])*T ))*mask[band])
                    plt.colorbar(fraction=.046,pad=.04)
                    plt.grid(color='w',ls='--')
                    plt.title('t_i - Tcorr')
                    plt.subplot(244)
                    plt.imshow((.5* p2[0][0] * ((2-p2[1][0])*T + p2[1][0] * np.cos(2*np.deg2rad(p2[2][0])) * Q + \
                                                p2[1][0] * np.sin(2*np.deg2rad(p2[2][0])) * U))*mask[band])
                    plt.title("Best Fit model")
                    plt.colorbar(fraction=.046,pad=.04)
                    plt.grid(color='w',ls='--')


                    plt.subplot(246)
                    plt.imshow((.5* p2[0][0] * ( p2[1][0] * np.cos(2*np.deg2rad(p2[2][0])) * Q )*mask[band]))
                    plt.title("Best Fit Q: %.2f"%(p2[1][0]* np.cos(2*np.deg2rad(p2[2][0]))))
                    plt.colorbar(fraction=.046,pad=.04)
                    plt.grid(color='w',ls='--')

                    plt.subplot(247)
                    plt.imshow((.5* p2[0][0] * ( p2[1][0] * np.sin(2*np.deg2rad(p2[2][0])) * U )*mask[band]))
                    plt.title("Best Fit U: %.2f"%(p2[1][0]* np.sin(2*np.deg2rad(p2[2][0]))))
                    plt.colorbar(fraction=.046,pad=.04)
                    plt.grid(color='w',ls='--')

                    plt.subplot(248)
                    residuals = (ti-.5* p2[0][0] * ((2-p2[1][0])*T + p2[1][0] * np.cos(2*np.deg2rad(p2[2][0])) * Q +\
                                                    p2[1][0] * np.sin(2*np.deg2rad(p2[2][0])) * U))*mask[band]

                    plt.imshow(residuals)
                    plt.title("t_i - Best-Fit")
                    plt.grid(color='w',ls='--')
                    plt.colorbar(fraction=.046,pad=.04)
                    plt.savefig('/home/axf295/2019/tmpPlots/%s_%sGHz_%sdeg.png'%(wafer,band,nomang))
                    '''
                    
                    i=0
                    #print('0th fit: ',p2)
                    while True:
                        res = CF.difference(p2[2][0],nomang,deg=True)
                        res0 = np.copy(res)
                        newnomang = nomang
                        while abs(res)>=30.:
                            i+=1
                            #print('firstfit',p2[2][0],nomang,res)
                            #print('%s may be mis-mapped. Rerunning fit_polarization_params with shifted nomang'%bolo)
                            last_res = np.copy(res)
                            p2=fit_polarization_params(nomang,band,ti,T,Q,U,mask[band],thetamean=p2[2][0],
                                                       dtheta=dtheta, gmean=p2[0][0],dg=dg,pmean=p2[1][0],
                                                       dp=dp,VARIABLE = variable)
                            res = CF.difference(p2[2][0],nomang,deg=True)
                            #print('%s iteration: '%i,p2)
                            #print('res:',res)
                            #print('last_res:',last_res)
                            if abs(last_res)<=abs(res):
                                res = last_res
                                break
        
                        if dtheta>1.:
                            i+=1
                            dtheta/=2.
                            dg /= 2.
                            dp /= 2.
                            p2=fit_polarization_params(nomang,band,ti,T,Q,U,mask[band],thetamean=p2[2][0],
                                                       dtheta=dtheta, gmean=p2[0][0],dg=dg,pmean=p2[1][0],
                                                       dp=dp,VARIABLE = variable)
                            res = CF.difference(p2[2][0],nomang,deg=True)
                            #print('%s iteration: '%i,p2)
                        else:
                            break
                    #print('best fit: ',p2)
        
                    
                    
                    ## rotate raw phi to be in [0,180]
                    ## Need 90 degree shift after changes to HWM on 5/30/2019
                    ## included fit to optimize rotation angle
                    
                    p2[2][0] -= rot_ang
                    if p2[2][0]<0:
                        p2[2][0]+=180.
                    if p2[2][0]>180.:
                        p2[2][0]-=180.
                       
                    phi = p2[2][0]
                    #print(p2)
                    ## For debugging; you can plot the different plots in the fitting process.
                    if PLOT:
                        plt.figure(figsize=(35,5))
                        plt.subplot(151)
                        plt.imshow(ti*mask[str(band)])
                        plt.title('%s, %sGHz, %s deg t_i'%(wafer,band,ang))
                        plt.colorbar()
                        plt.grid(color='w',ls='--')
                        plt.subplot(152)
                        plt.imshow(T*mask[str(band)])
                        plt.colorbar()
                        plt.grid(color='w',ls='--')
                        plt.title('T')
                        plt.subplot(153)
                        plt.imshow((ti-T)*mask[str(band)])
                        plt.colorbar()
                        plt.grid(color='w',ls='--')
                        plt.title('t_i-T')
                        plt.subplot(154)
                        plt.imshow((.5* p2[0][0] * ((2-p2[1][0])*T + p2[1][0] * np.cos(2*np.deg2rad(p2[2][0])) * Q + p2[1][0] * np.sin(2*np.deg2rad(p2[2][0])) * U))*mask[str(band)])
                        plt.title("Best Fit t_i")
                        plt.colorbar()
                        plt.grid(color='w',ls='--')
                        plt.subplot(155)
                        residuals = ((.5* p2[0][0] * ((2-p2[1][0])*T + p2[1][0] * np.cos(2*np.deg2rad(p2[2][0])) * Q + p2[1][0] * np.sin(2*np.deg2rad(p2[2][0])) * U))-ti)*mask[str(band)]
                        residuals = residuals[np.nonzero(residuals)]
                        dx = np.shape(ti)[0]-20
                        var = np.nanvar(ti[dx:,:20])
                        plt.hist(np.ndarray.flatten(residuals)/np.sqrt(var),bins=np.arange(-5,5,.1))
                        x = np.arange(-5,5,.01)
                        y = 25.*np.exp(-x**2/2.)
                        plt.plot(x,y,color='r')
                        plt.title("Fit Residual / var")
                        plt.grid(color='w',ls='--')
                        plt.savefig(plot_loc+'%s_%sGHz_%sdeg_FittingPlot.png'%(wafer,band,ang))
                        plt.close('all')

                    nomang =  ang
                    residual_angle = CF.difference(phi,nomang,deg=True)
                    
                    
                    log.write('Best-fit Phi: %.2f\n'%phi)
                    log.write('Raw phi: %.2f +- %.2f\n'%(p2[2][0],p2[2][1]))
                    log.write('Best-fit Gain: %.2f +- %.2f\n'%(p2[0][0],p2[0][1]))
                    log.write('Best-fit PolGain: %.2f +- %.2f\n'%(p2[1][0],p2[1][1]))

                    fitdata[wafer][ang][band]['nbolos'] = numbolos[wafer][ang][band]
                    fitdata[wafer][ang][band]['G']      = p2[0]
                    fitdata[wafer][ang][band]['p']      = p2[1]
                    fitdata[wafer][ang][band]['rawPhi'] = p2[2]
                    fitdata[wafer][ang][band]['nom']    = nomang
                    fitdata[wafer][ang][band]['rotPhi'] = phi
                    fitdata[wafer][ang][band]['res']    = residual_angle
                    

    return fitdata
