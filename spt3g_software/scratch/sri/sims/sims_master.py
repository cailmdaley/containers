import numpy as np, scipy as sc, os, sys, time, glob
from scipy.io import readsav
from scipy import optimize
import healpy as H
from pylab import *

class sims_master():

	def __init__(self):

		"""
		set all default variables for sims here
		"""

		self.spt3g_software_PATH_envvar = 'spt3g_software_PATH'

		#main folder
		#change
		##self.parent_folder = '/Users/sraghunathan/Research/SPTPol/analysis/2018_10/spt3g_sims'

		
		##self.data_folder = '/Users/sraghunathan/Research/SPTPol/analysis/2018_10/spt3g_sims/data/' #'/spt/user/sri/data/'
		##self.data_folder = '/home/sri/analysis/2018_10/spt3g_sims/data/'
		try:			
			self.data_folder = '%s/scratch/sri/data/' %(os.environ[self.spt3g_software_PATH_envvar])
		except:
			self.data_folder = '/home/sri/spt3g_software/scratch/sri/data/'

		#units
		self.Tcmb = 2.725 #CMB temeprature in kelvin
		self.K_to_uK = 1e6 #kelvin to uK

		#degrees/radians
		self.degrees2radians = np.pi / 180.
		self.arcmins2radians = self.degrees2radians / 60.

		#CAMB details
		##self.camb_path = '/Users/sraghunathan/Research/SPTPol/analysis/2018_10/spt3g_sims/tools/camb' ##'/spt/user/sri/tools/camb/'
		##self.camb_path = '/home/sri/analysis/2018_10/spt3g_sims/tools/camb/'
		##self.camb_obj_module = 'camb'
		self.input_cosmology_for_camb = 'base_plikHM_TTTEEE_lowl_lowE_lensing'
		try:			
			self.camb_path = '%s/simulations/camb' %(os.environ[self.spt3g_software_PATH_envvar])
		except:			
			self.camb_path = '/home/sri/spt3g_software/simulations/camb'
			print('\n\t\t%s path not set as env. variable. reading CAMb files from %s') %(self.spt3g_software_PATH_envvar, self.camb_path)

		self.reqd_camb_output_field = 'scalCls' #only unlensed CAMB power spectra now

		#verbose
		self.verbose = 0

		#beams
		###self.gaussian_beam_dic = {90: 1.7, 150: 1.2, 220: 1.0}

		#map stuff
		self.nside = 256 #healpix map resolution
		self.lmax = 768 #max \ell for healpix operations 
		self.use_mask = 1 #apply 3G mask
		self.do_pol = 1 #do T/Q/U
		self.spt3g_mask_file = '/home/sri/3g_stuffs/MASK.fits'

		#foreground stuff
		self.SZ_pol_fraction = 0.04 #pol. fraction for SZ
		self.DG_pol_fraction = 0.02 #pol.fraction of dusty galaxies
		self.RG_pol_fraction = 0.04 #pol. fraction of radio galaxies
		self.george_fg_dic_file = '%s//foreground/george_plot_bestfit_line.sav' %(self.data_folder)
		self.de_zotti_number_counts_file = '%s//foreground/counts_150GHz_de_Zotti_radio.res' %(self.data_folder)
		self.dust_spec_index_radio_george_spt15 = -0.9 #to scale dust from 150 GHz
		self.bethermin_number_counts_file_prefix = '%s//foreground/Bethermin_model_counts_' %(self.data_folder)
		self.flux_cut_limit = 50e-3 #mJy based on G15 (is there a max cut - not quite sure.)
		self.min_flux_limit = 6.4e-3 #mJy based on G15 (i think this is the 5\sigma limit for G15)


		self.separator = '\n############################'

	def is_seq(self,o):
		return hasattr(o, '__len__')

	def fn_get_units(self, G3Units = 1):

		"""
		get map units here
		"""
		if G3Units:
			import spt3g.core
			units = core.G3Units.uK  # ensure maps are stored with correct unit information
		else:
			units = self.K_to_uK

	def fn_check_input_params(self, param, paramval): #TBC

		"""
		make sure the param values make sense		
		"""

		to_check_dic = {'which_cosmology': ['planck_2015_conservative']}

		if param in to_check_dic: 
			if paramval not in to_check_dic[param]: 
				return 0
			else: 
				return 1

		return 1

	def fn_get_param_dict_from_paramfile(self, paramfile, copy_to_class = 0, validate_params = 0):

		"""
		this function reads the paramfile, validates it if required, and stores/returns it as a dictionary: self.param_dict
		"""

		params, paramvals = np.genfromtxt(paramfile, delimiter = '=', unpack = 1, dtype = None, autostrip = True)
		param_dict = {}
		for p,pval in zip(params,paramvals):

			try:
				pval = float(pval)
			except:
				pass

			if (1):#replace unallowd characters in paramname
				p = p.replace('(','').replace(')','')

			#check if the param value is valid #TBC
			if validate_params: 
				valid_param = self.fn_check_input_params(p,pval)
				if not valid_param: continue

			#push it into dic
			##from IPython import embed; embed()
			param_dict[p] = pval

		#save it as a class variable
		if copy_to_class: self.param_dict = param_dict

		#overwrite default variables in class with values from param_dict
		for keyname in param_dict:
			if isinstance(param_dict[keyname], str):
				cmd = 'self.%s = "%s"' %(keyname, param_dict[keyname])
			else:
				cmd = 'self.%s = %s' %(keyname, param_dict[keyname])
			exec(cmd)

		return param_dict

	def fn_sim_setup(self, param_dict):

		"""
		set up sim stuff
		step 1. use the CAMB params file to check if CAMB must be run
		step 2. read an existing power spectra or the currently generated power spectra
		"""

		'''#TBC: Not running CAMB - we will just read CAMB output file
		#run CAMB
		self.fn_run_camb(param_dict)
		'''

		#read CAMB output
		self.els, self.Cls = self.fn_read_camb(param_dict)

	def fn_file_exists(self, fname):
		return os.path.exists(fname)

	def fn_create_subdirs(self, opfolder):
		success = 0
		opfolder_split = opfolder.split('/')
		tmp = ''
		for fcnt, ff in enumerate( opfolder_split ):
			tmp = '%s/%s' %(tmp,ff);
			if tmp =='/': continue
			if fcnt == 0: tmp = tmp.strip('/')
			if not os.path.exists(tmp):
				success = os.system('mkdir %s' %(tmp))
				if success>0:
					print('\n\t\tAlert: %s could not be created' %(opfolder))
		return success

	def fn_get_gauss_beam(self, nu, lmax):
		beamfwhm = self.gaussian_beam_dic[nu]
		beamfwhm *= self.arcmins2radians
		
		return H.gauss_beam(beamfwhm, lmax)


	'''#TBC: commenting LOC that will run CAMB
	def fn_get_camb_power_spectra_fname(self, param_dict, camb_param_file):

		camb_dict = self.fn_get_param_dict_from_paramfile(camb_param_file)
		camb_output_file_searchstr = '%s*%s.dat' %(camb_dict['output_root'], param_dict['reqd_camb_output_field'])
		camb_output_file = glob.glob(camb_output_file_searchstr)

		if len(camb_output_file) == 0:
			return None
		else:
			return camb_output_file[0]

	def fn_get_camb_param_fname(self, input_cosmology_for_camb):

		camb_param_file_prefix = '%s/params' %(self.camb_path)
		camb_param_file_suffix = 'params_%s.ini' %(input_cosmology_for_camb)
		camb_param_file = '%s/%s' %(camb_param_file_prefix, camb_param_file_suffix)

		return camb_param_file

	def fn_run_camb(self, param_dict):

		"""
		this function determines if we must run CAMB for the specified input or just read an already existing CAMB power spectra
		"""

		camb_param_file = self.fn_get_camb_param_fname(param_dict['input_cosmology_for_camb'])
		print camb_param_file
		sys.exit()

		#check if the param file for CAMB exists
		if not self.fn_file_exists(camb_param_file): 
			logline = '\n\tAborting: %s file to execute CAMB is not found\n' %(camb_param_file)
			print(logline)
			sys.exit()

		#check if output already exists
		run_camb = 0
		camb_output_file = self.fn_get_camb_power_spectra_fname(param_dict, camb_param_file)
		if camb_output_file == None: run_camb = 1

		#if outputs already exist just use them
		if not run_camb: return

		#else run CAMB now
		#make sure output folder exists. else create them
		camb_dict = self.fn_get_param_dict_from_paramfile(camb_param_file)
		cambopfolder = '/'.join(camb_dict['output_root'].split('/')[0:-1])
		self.fn_create_subdirs(cambopfolder)
		camb_run_cmd = '%s/%s %s' %(self.camb_path, self.camb_obj_module, camb_param_file)
		print('\trunning CAMB now')
		os.system(camb_run_cmd)

		#ensure output file is created
		camb_output_file = self.fn_get_camb_power_spectra_fname(param_dict, camb_param_file)		
		if not self.fn_file_exists(camb_output_file):
			logline = '\n\tAlert: %s output file from CAMB is not found. Something is wrong with CAMB. Check\n' %(camb_output_file)
			print(logline)
			sys.exit()

		return

	def fn_read_camb(self, param_dict):

		"""
		this function read the CAMB output power spectra file to produce sims
		"""

		camb_param_file = self.fn_get_camb_param_fname(param_dict['input_cosmology_for_camb'])
		camb_output_file = self.fn_get_camb_power_spectra_fname(param_dict, camb_param_file)		
		if not self.fn_file_exists(camb_param_file): 
			logline = '\n\tAborting: %s output file from CAMB is not found. Run CAMB to get it.\n' %(camb_param_file)
			print(logline)
			self.fn_run_camb(param_dict)

		#get the power spectra now
		els, Cls = self.fn_get_camb_power_spectra(camb_output_file, camb_param_file)

		return els, Cls
	'''#TBC: commenting LOC that will run CAMB

	def fn_get_camb_power_spectra_fname(self, param_dict, camb_param_file):

		"""
		reads the CAMB ini to get the output_root for CAMB Cls file
		"""

		camb_dict = self.fn_get_param_dict_from_paramfile(camb_param_file)
		camb_output_file_searchstr = '%s/%s*%s.dat' %(self.camb_path, camb_dict['output_root'], param_dict['reqd_camb_output_field'])
		camb_output_file = glob.glob(camb_output_file_searchstr)

		try:
			return camb_output_file[0]
		except:
			return None

	def fn_get_camb_param_fname(self, input_cosmology_for_camb):

		"""
		gets the CAMB ini file
		"""

		##camb_param_file = '%s/%s_for3g.ini' %(self.camb_path, input_cosmology_for_camb)
		##return camb_param_file

		camb_param_file_arr = glob.glob('%s/%s*.ini' %(self.camb_path, input_cosmology_for_camb))
		try:
			return camb_param_file_arr[0]
		except:
			return None		

	def fn_read_camb(self, param_dict):

		"""
		this function read the CAMB output power spectra file to produce sims
		"""

		#first get the CAMB params file
		camb_param_file = self.fn_get_camb_param_fname(param_dict['input_cosmology_for_camb'])

		#now get the CAMB Cls file
		camb_output_file = self.fn_get_camb_power_spectra_fname(param_dict, camb_param_file)

		if not self.fn_file_exists(camb_output_file): 
			logline = '\n\tAborting: %s output file from CAMB is not found. Run CAMB to get it.\n' %(camb_param_file)
			print(logline)
			self.fn_run_camb(param_dict)

		#get the power spectra now
		els, Cls = self.fn_get_camb_power_spectra(camb_output_file, camb_param_file)

		return els, Cls

	def fn_get_camb_power_spectra(self, camb_output_file, camb_param_file = None):

		"""
		this function outputs the CAMB output power spectra
		"""

		Dls = np.loadtxt(camb_output_file, usecols=[1,2,3,4])
		els = np.loadtxt(camb_output_file, usecols=[0])

		if camb_param_file <> None: #Check if Tcmb is included for not from CAMB ini file
			#take the Tcmb factor into account here
			camb_dict = self.fn_get_param_dict_from_paramfile(camb_param_file)
			if float(camb_dict['CMB_outputscale']) == 1.:
				Dls = camb_dict['temp_cmb']**2. * Dls

		Dls_fac = ( els[:,None] * (els[:,None]  + 1) ) / (2 * np.pi )

		Cls = Dls / Dls_fac

		Cls = Cls.T

		return els, Cls

	def fn_synfast(self, Cls, nside, lmax = None, fwhm = 0., sigma = None):

		"""
		make sims using synfast
		"""

		return H.synfast(Cls, nside, lmax = lmax, fwhm = fwhm, sigma = sigma, verbose = self.quiet)

	def fn_get_set_random_seed(self, sim_index, nu, set_seed = 1):

		"""
		get and set random seed
		"""

		random_seed = int( (sim_index + 1 ) * nu )

		if set_seed: self.fn_set_random_seed(random_seed)

		return random_seed

	def fn_set_random_seed(self, random_seed):

		"""
		set a given random seed
		"""

		return np.random.seed(random_seed)

	def fn_Jy_K(self, nu, Jy_K = 1, K_Jy = 0):

		"""
		Jy to K conversion and vice versa
		"""

		nu *= 1e9
		h, k, c, temp=6.626e-34,1.38e-23, 3e8, 2.725
		x=h*nu/(k*temp)

		dbnu = 2.*k*(nu**2/c**2)*(x**2*np.exp(x))/(np.exp(x)-1.)**2.0
		conv=1e26*dbnu #Jy to K

		if Jy_K:
			return 1./conv
		else:
			return conv


	def fn_get_foreground_power_George_SPT_2015(self, which_fg, nu1 = 150, nu2 = None, ign_DG_RG = 1, check_fg_labels = 1):

		"""
		foreground powers from George_SPT_2015 results
		"""

		which_fg = which_fg.replace('_Gaussian','')
		possible_fgs = ['all', 'tSZ', 'kSZ','DG-Cl','DG-Po','RG','tSZ-CIB','Total','CMB']
		#possible_fgs = ['all','DG-Cl','DG-Po']
		if check_fg_labels: assert which_fg in possible_fgs

		try:		
			self.george_fg_dic
		except:
			self.george_fg_dic = readsav(self.george_fg_dic_file)

		if nu2 == None: nu2 = nu1

		if nu1 == 90: nu1 = 95
		if nu2 == 90: nu2 = 95
		
		fg_nus = np.asarray( [(95,95), (95,150), (95,220), (150,150), (150,220), (220,220)] )

		fg_labels = self.george_fg_dic['ml_dl_labels']
		fg_nu_ind = np.where((fg_nus[:,0] == nu1) & (fg_nus[:,1] == nu2))[0][0]

		fg_els = self.george_fg_dic['ml_l']

		if which_fg <> 'all':
			fg_lab_ind = np.where(fg_labels == which_fg)[0]
			fg_Dls = self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0]
		else:
			fg_Dls_all = self.george_fg_dic['ml_dls'][fg_nu_ind]
			for fgcnt, fg in enumerate(possible_fgs):

				########################
				if fg == 'all'  or fg == 'tSZ-CIB': continue

				fg_lab_ind = np.where(fg_labels == fg)[0]

				if ign_DG_RG and (fg == 'DG-Cl' or fg == 'DG-Po' or fg == 'RG'): continue
				if fg == 'Total' or fg == 'CMB': continue
				########################
				try:
					fg_Dls += self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0]
				except:
					fg_Dls = self.george_fg_dic['ml_dls'][fg_nu_ind][fg_lab_ind][0]

		fg_Cls = ( fg_Dls * 2 * np.pi ) / ( fg_els * (fg_els + 1) ) * self.K_to_uK**-2.

		return fg_Cls

	def fn_get_pol_fraction(self, param_dict, which_fg):

		"""
		set up polarisation fraction		
		"""

		try:
			if which_fg.find('SZ')>-1:
				return param_dict['SZ_pol_fraction']
			elif which_fg.find('DG')>-1:
				return param_dict['DG_pol_fraction']
			elif which_fg.find('RG')>-1:
				return param_dict['RG_pol_fraction']
		except:
			return self.pol_fraction

	def fn_set_up_foreground_sims(self, param_dict, nu):

		"""
		get the required foreground power spectra
		"""

		FG_Cls_dic = {}
		#for Gaussian realisation of all foregrounds
		if param_dict['all_Gaussian']:
			print('\t\tgetting all FG power spectra from George et al. 2015 measurements')
			all_Gaussian_FG_Cls = self.fn_get_foreground_power_George_SPT_2015(which_fg, nu1 = nu, nu2 = nu)
			FG_Cls_dic['all_Gaussian_FG_Cls'] = all_Gaussian_FG_Cls

		#for Gaussian SZ
		if param_dict['make_SZ']:
			print('\t\tgetting tSZ power spectra from George et al. 2015 measurements')
			tSZ_Cls = self.fn_get_foreground_power_George_SPT_2015('tSZ', nu1 = nu, nu2 = nu)
			FG_Cls_dic['tSZ'] = tSZ_Cls

			print('\t\tgetting kSZ power spectra from George et al. 2015 measurements')
			kSZ_Cls = self.fn_get_foreground_power_George_SPT_2015('kSZ', nu1 = nu, nu2 = nu)
			FG_Cls_dic['kSZ'] = kSZ_Cls

		#for Gaussian radio / dust
		if (1):##param_dict['Gaussian_psources']:

			#radio Poisson component
			if param_dict['make_RG']:
				print('\t\tgetting RG power spectra from George et al. 2015 measurements')
				RG_Cls = self.fn_get_foreground_power_George_SPT_2015('RG', nu1 = nu, nu2 = nu)
				FG_Cls_dic['RG'] = RG_Cls

			#dusty galaxies Poisson component
			if param_dict['make_DG']:
				print('\t\tgetting DG Poisson term power spectra from George et al. 2015 measurements')
				DGPo_Cls = self.fn_get_foreground_power_George_SPT_2015('DG-Po', nu1 = nu, nu2 = nu)
				FG_Cls_dic['DG'] = DGPo_Cls
				if (1):
					DGclus_Cls = self.fn_get_foreground_power_George_SPT_2015('DG-Cl', nu1 = nu, nu2 = nu)
					FG_Cls_dic['DG'] += DGclus_Cls
		
		#dusty galaxies clustering term. get this even if we do Poisson realisations to include the clustering term later
		if param_dict['make_DG_clus']:
			print('\t\tgetting DG clustering term power spectra from George et al. 2015 measurements')
			DGclus_Cls = self.fn_get_foreground_power_George_SPT_2015('DG-Cl', nu1 = nu, nu2 = nu)
			FG_Cls_dic['DGclus'] = DGclus_Cls

		self.FG_Cls_dic = FG_Cls_dic

		return FG_Cls_dic

	def fn_get_Poisson_source_counts(self, which_pop, nu = 150, do_pol = 1, below_flux_limit = 1, assign_fluxes = 1, dust_spec_index_radio = None):


		"""
		get dusty or radio source counts for Poisson realisaitons
		flux_limit is the limit for point source masking [in Jy]
		min_flux_limit also in Jy: set to 1mJy
		We will simulated all sorces between 1 and 5 mJy. Anything above that should be masked in the maps
		"""

		assert nu == 95 or nu == 150 or nu == 220
		assert which_pop == 'RG' or which_pop == 'DG' #radio or dusty galaxies

		"""
		#first we will get flux range [s] and dN/ds for dustry or radio sources
		"""

		if which_pop == 'RG':

			#pending - handle different models other than just deZotti
			de_zotti_number_counts = np.loadtxt(self.de_zotti_number_counts_file, skiprows = 12)
			radio_log_s,radio_x,radio_y,radio_logs25n = np.asarray( zip(*de_zotti_number_counts) )


			# first get the number of sources in each flux range
			s =  10**radio_log_s

			#scale to other frequencies from 150 GHz
			if nu <> 150:
				if dust_spec_index_radio == None: 
					dust_spec_index_radio = self.dust_spec_index_radio_george_spt15 #goerge et al. 2015
				s = (150./nu)**dust_spec_index_radio * s

			s25n = 10**radio_logs25n
			dnds = s25n / s**2.5

		elif which_pop == 'DG':

			#pending - handle different models other than just Bethermin
			if nu == 95:
				nu_bethermin = 100
			elif nu == 150:
				nu_bethermin = 143
			elif nu == 220:
				nu_bethermin = 217

			bethermin_number_counts_file = '%s%sGHz.txt' %(self.bethermin_number_counts_file_prefix, nu_bethermin)
			bethermin_number_counts  = np.loadtxt(bethermin_number_counts_file,usecols = [0,1])

			s,dnds = np.asarray( zip(*bethermin_number_counts) )

		lns = np.log(s)
		dlns =  (lns[1]-lns[0])
		dndlns = dnds*s
		nsources=(dndlns)*dlns #number of sources obtained

		#pick sources below or above flux cut
		if below_flux_limit: 
			sinds = np.where((s<self.flux_cut_limit))[0]
		else:
			#pick sources above a certain flux
			sinds = np.where((s>self.flux_cut_limit))[0]

		s, nsources = s[sinds], nsources[sinds]

		if (1):#just makes the program run faster by ignoring extremely faint sources below 1e-7 J
			sinds = np.where((s>self.min_flux_limit))[0]
			s, nsources = s[sinds], nsources[sinds]

		psource_counts = [s, nsources] #[flux bin, no of sources in each bin]

		return psource_counts 

	def fn_make_Poisson_source_realisation(self, which_pop, param_dict, nu, sim_index, noofsims, psource_counts = None, random_seed = -1, nside = None, lmax = None, weights = None, show_plot = 0):

		if not self.is_seq(psource_counts):
			psource_counts = self.fn_get_Poisson_source_counts(which_pop, nu = 150)

		#nside
		if nside == None: nside = int( param_dict['nside'] )

		#lmax
		if lmax == None: lmax = int( param_dict['lmax'] )

		s, nsources = psource_counts

		dx = H.nside2resol(nside)
		npixels = H.nside2npix(nside)

		pix_area = dx ** 2. #Sr
		lam = nsources * pix_area#; print np.sum(lam)

		#weights - any pixel based weights? DG clustering or tSZ x CIB correlation?
		if not self.is_seq(weights): weights = np.ones(npixels)
		weights = weights/np.sum(weights)

		#do_pol: just T or TQU
		try:
			do_pol = param_dict['do_pol']
		except:
			do_pol = 0

		#random_seed
		if random_seed <> -1: 
			self.fn_set_random_seed(random_seed)
		else:
			if noofsims == 1: #then set random seed based on a specific sim_index and nu
				random_seed = self.fn_get_set_random_seed(sim_index, nu )

		#Jy to K
		Jy_to_K_conv = self.fn_Jy_K(nu)

		if param_dict['use_mask']: #apply 3G mask
			MASK = H.ud_grade( H.read_map('MASK.fits', verbose = self.verbose), nside)

		if do_pol: 
			pol_fraction = self.fn_get_pol_fraction(param_dict, which_pop)

		SIMS = []
		for simcnt in range(noofsims): #loop over simulations now

			###TBC: currently i am stupidly doing this for each flux bin. Need to be fixed using the below code where i can do this for all bins at once and pick flux values from the dN/dS distribution

			logds = np.log10(s[1]) - np.log10(s[0])
			logs_range = np.asarray([(np.log10(val), np.log10(val) + logds) for val in s])
			s_range = 10**logs_range


			POIS_SIMS_IN_BINS= []
			if do_pol: #create Q, U maps too
				POIS_SIMS_Q_IN_BINS = []
				POIS_SIMS_U_IN_BINS = []

			for scnt, (s1,s2) in enumerate(s_range): #create a source map for each flux bin

				###print s1, s2

				lam_vec = np.zeros(npixels) + lam[scnt]
				lam_vec = (np.sum(lam_vec) * weights)##.flatten()

				#createing number counts map
				CURR_BIN = np.asarray( map(lambda x: np.random.poisson(x, 1) * 1., lam_vec) ).T[0]
				if do_pol: #create Q, U maps too
					CURR_BIN_Q = np.zeros( len(CURR_BIN) )
					CURR_BIN_U = np.zeros( len(CURR_BIN) )
				##print CURR_BIN, np.sum(CURR_BIN)

				if show_plot:
					H.mollview( CURR_BIN, title = 'Number counts' ); show()

				if param_dict['use_mask']: #apply 3G mask, if need be
					CURR_BIN *= MASK
					if do_pol:
						CURR_BIN_Q *= MASK
						CURR_BIN_U *= MASK

					if show_plot:
						H.mollview(CURR_SIM, title = 'number counts after SPT-3G masking');show();sys.exit()

				## now assign fluxes
				NON_ZERO_INDS = np.where(CURR_BIN>0)[0]
				for nn in NON_ZERO_INDS:

					#we can approximate the flux value to folow a uniform distribution in each narrow flux bin
					fluxes = np.random.uniform(s1,s2,size=int(CURR_BIN[nn])) #in Jy now
					CURR_BIN[nn] = np.sum(fluxes) * Jy_to_K_conv / pix_area #in K now

					if do_pol:
						P = fluxes * pol_fraction
						pol_angle = np.random.uniform(0.,360.,len(fluxes)) #a random pol angle for each source
						Q = P * np.cos(np.radians(2*pol_angle))
						U = P * np.sin(np.radians(2*pol_angle))

						CURR_BIN_Q[nn] = np.sum(Q) * Jy_to_K_conv / pix_area #in K now
						CURR_BIN_U[nn] = np.sum(U) * Jy_to_K_conv / pix_area #in K now

				if show_plot:
					H.mollview( CURR_BIN_Q, title = 'Poisson source power' ); show()

				POIS_SIMS_IN_BINS.append(CURR_BIN)
				if do_pol:
					POIS_SIMS_Q_IN_BINS.append( CURR_BIN_Q )
					POIS_SIMS_U_IN_BINS.append( CURR_BIN_U )

			#adding the flux values obtained from all bins
			CURR_SIM = np.sum( POIS_SIMS_IN_BINS, axis = 0 )
			if do_pol:
				CURR_SIM_Q = np.sum( POIS_SIMS_Q_IN_BINS, axis = 0 )
				CURR_SIM_U = np.sum( POIS_SIMS_U_IN_BINS, axis = 0 )


			if show_plot:
				##from IPython import embed; embed()
				H.mollview( CURR_SIM[0] ); show();sys.exit()

			if do_pol:
				SIMS.append( [CURR_SIM, CURR_SIM_Q, CURR_SIM_U] )
			else:
				SIMS.append( CURR_SIM )

		return np.asarray(SIMS)

		"""
		if (1):

			for simcnt in range(noofsims): #loop over simulations now
				lam_vec = np.zeros(npixels) + np.sum(lam)
				CURR_SIM = np.asarray( map(lambda x: np.random.poisson(x, 1)[0], lam_vec))
				CURR_SIM = CURR_SIM * 1.
				##H.mollview(CURR_SIM, title = 'number counts before masking');show();sys.exit()

				if use_mask or param_dict['use_mask']: #apply 3G mask
					MASK = H.ud_grade( H.read_map('MASK.fits', verbose = self.verbose), nside)
					CURR_SIM = CURR_SIM * MASK
					###H.mollview(CURR_SIM, title = 'number counts after SPT-3G masking');show()

				for n in range(npixels):
					if int( CURR_SIM[n] ) == 0: continue
					if n%5000 == 0: print n, npixels
					fluxes = self.fn_random_sampler(s, lam, howmanysamples = int( CURR_SIM[n] ), burn_in = 0)
					from IPython import embed; embed()

					min_edge=min(np.log10(s)); max_edge=max(np.log10(s))
					binbin=np.logspace(min_edge, max_edge, len(s))
					ax = subplot(111, yscale = 'log', xscale = 'log')
					hist(fluxes, bins = s);plot(s, lam, color = 'orange')
					show()

					##print np.sum(fluxes)* Jy_to_K_conv / pix_area
					CURR_SIM[n] = np.sum(fluxes) * Jy_to_K_conv / pix_area #in K now
					
				from IPython import embed; embed()
				H.mollview(CURR_SIM * 1e6, title = 'Poisson sources - temperature map');show()


				if (1):				
					clf(); ax = subplot(111, yscale = 'log')
					hist(fluxes, bins = s, normed = 0); plot(s, lam, 'r')
					show();sys.exit()

				print CURR_SIM;sys.exit()
		"""


	def fn_make_Gaussian_simulated_skies(self, param_dict, sim_index, nu = 150, cmb_sky = 1, fg_sky = 0, noofsims = 1, random_seed = -1, Cls = None, nside = None, lmax = None, do_pol = 1):

		"""
		creates Gaussian CMB or foreground sims using healpix routines
		"""

		if not self.is_seq(Cls): 
			if cmb_sky: 
				print('\tgetting required CMB power spectra from CAMB')
				Cls = self.Cls #get Cls first from CAMB
			elif fg_sky:
				print('\tgetting required foreground power spectra')
				FG_Cls_dic = self.fn_set_up_foreground_sims(param_dict, nu)


		#nside
		if nside == None: nside = int( param_dict['nside'] )

		#lmax
		if lmax == None: lmax = int( param_dict['lmax'] )

		#do_pol: just T or TQU
		try:
			do_pol = param_dict['do_pol']
		except:
			do_pol = 0

		#random_seed
		if random_seed <> -1: 
			self.fn_set_random_seed(random_seed)
		else:
			if noofsims == 1: #then set random seed based on a specific sim_index and nu
				random_seed = self.fn_get_set_random_seed(sim_index, nu )

		'''
		#get beam
		if param_dict['exp_beam'] == 0: #Healix Gaussian beam
			Bl = self.fn_get_gauss_beam(nu, lmax)
			###plot(Bl);show();sys.exit()
		'''

		#create simulated CMB/foreground skies
		"""
		dimension of SIM_MAPS will be ( noofsims x teb_len x npixels )
		noofsims = total number of simulations
		teb_len = 3 #if do_pol = 1 #else teb_len = 1
		npixels = map_pixels based on healpix nside
		"""
		if int(param_dict['use_mask']): #apply 3G mask
			MASK = H.ud_grade( H.read_map(self.spt3g_mask_file, verbose = self.verbose), nside)

		if cmb_sky:
			SIMS = np.asarray( [H.synfast(Cls, nside, lmax = lmax, verbose = self.verbose, pol = do_pol) for n in range( noofsims )] )
			if param_dict['use_mask']: #apply 3G mask
				SIMS *= MASK
			return SIMS

		if fg_sky:

			FG_SIM_MAPS_dic = {}
			for keyname in FG_Cls_dic.keys():
				if keyname == 'DGclus': continue
				curr_FG_Cls = FG_Cls_dic[keyname]

				SIMS = np.asarray( [H.synfast(curr_FG_Cls, nside, lmax = lmax, verbose = self.verbose, pol = do_pol) for n in range( noofsims )] )
				if do_pol: 
					pol_fraction = self.fn_get_pol_fraction(param_dict, keyname)
					SIMS = np.asarray( map(lambda x: [x, x*pol_fraction, x*pol_fraction], SIMS) )

				if param_dict['use_mask']: #apply 3G mask
					SIMS *= MASK

				FG_SIM_MAPS_dic[keyname] = SIMS

			return FG_SIM_MAPS_dic

	def fn_dump_files(self, SIMS, parentfolder, subfolder, sim_pref = 'sim'):
		"""
		dumps the simulations in the respective folder
		#TBC: not yet complete - we will finish this based on what needs to be stored
		"""
		opfoldersuccess = self.fn_create_subdirs(parentfolder) #create folders for outputs - currently no output is stored
		if opfoldersuccess == 0: #now dump the files (not completed yet. depends on what we want to store)
			simsfoldername = '%s/%s' %(parentfolder, subfolder)
			success = self.fn_create_subdirs(simsfoldername)
			if success == 0: #store sims now
				noofsims = len(SIMS)
				for n in range(noofsims):
					cmbfname = '%s/%s_%s.fits' %(simsfoldername, sim_pref, n)
					print(cmbfname)
