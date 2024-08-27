import numpy as np

class Spectrum(object):
    """
    Spectrum class to store individual FAS object (containing also info on 
    the event originating it).
    The class features are:
      freq          <--> [array, N] frequency points
      amp           <--> [array, N] FAS amplitude values
      sm_freq       <--> [array, M] undersampled frequency points
      sm_amp        <--> [array, M] smoothed FAS amplitude values
      fbool         <--> [array, M] boolean values to flag datapoints used
                          for calculations (on smoothed spectra)
      delta         <--> [scalar] source-station angular distance
      orid          <--> [scalar] event identifier
      sta           <--> [scalar] station code 
      chan          <--> [scalar] channel code
      hypdist       <--> [scalar] hypocentral distance in km
      ml            <--> [scalar] local magnitude of the event
      mag           <--> [scalar] moment magnitude of the event
      freq_lims     <--> [array, 2] lower and upper frequencies used for 
                          calculations
    """

    def __init__(self, freq=np.array([]), amp=np.array([]), 
                 sm_freq=np.array([]), sm_amp=np.array([]),
                 fbool=np.array([]), delta=[], orid=None,
                 sta=None, chan=None, hypdist=None, ml=None, mag=None,
                 freq_lims=np.array([0., 0.])):

        self.freq = freq
        self.amp = amp
        self.sm_freq = sm_freq
        self.sm_amp = sm_amp
        self.fbool = fbool
        self.delta = delta
        self.orid = orid
        self.sta = sta
        self.chan = chan
        self.hypdist = hypdist
        self.ml = ml
        self.mag = mag
        self.freq_lims = freq_lims


    def __str__(self):
        string = u'[<Spectrum> orid:%s sta:%s chan:%s hypdist:%s ml:%s ' \
                 u'f1:%s f2:%s delta:%s \n sm_freq:%s \n'\
                 u'sm_amp:%s \n fbool:%s \n' % (self.orid, self.sta, self.chan,
                                     self.hypdist, self.ml, self.freq_lims[0],
                                     self.freq_lims[1], self.delta, 
                                     self.sm_freq, self.sm_amp, self.fbool, )
        return string



class StaticParams(object):
    """
    Static parameters class.
    The class features are:
      N_ev          <--> [scalar] number of events
      N_sta         <--> [scalar] number of stations
      N_freq        <--> [scalar] number of frequency points (in each spectrum)
      M             <--> [array, N_ev*N_sta*N_freq] boolean values to flag
                          datapoints used for calculations
      F             <--> [array, N_ev*N_sta*N_freq] frequency points
      R             <--> [array, N_ev*N_sta*N_freq] hypocentral distances
      data          <--> [array, N_ev*N_sta*N_freq] FAS data
      weights       <--> [array, N_ev*N_sta*N_freq] optional weights applied
                          to data
      scale         <--> [scalar] optional scaling for cost function
                          calculation and for forward modelling
      orids         <--> [array, N_ev] list of event origin IDs
      stas          <--> [array, N_sta] list of stations
      freqs         <--> [array, N_freqs] list of frequency points in each
                          spectrum
    """    

    def __init__(self, M=None, F=None, R=None, N_ev=None, N_sta=None, N_freq=None, data=None, weights=None, scale=None, orids=None, stas=None, freqs=None):

        self.M = M
        self.F = F
        self.R = R
        self.N_ev = N_ev
        self.N_sta = N_sta
        self.N_freq = N_freq
        self.data = data
        self.weights = weights
        self.scale = scale
        self.orids = orids
        self.stas = stas
        self.freqs = freqs

    def set_M(self, M):
        self.M = M

    def set_F(self, F):
        self.F = F
        
    def set_R(self, R):
        self.R = R

    def set_N_ev(self, N_ev):
        self.N_ev = N_ev

    def set_N_sta(self, N_sta):
        self.N_sta = N_sta

    def set_N_freq(self, N_freq):
        self.N_freq = N_freq

    def set_data(self, data):
        self.data = data

    def set_scale(self, scale):
        self.scale = scale

    def set_weights(self, weights):
        self.weights = weights

    def set_orids(self, orids):
        self.orids = orids

    def set_stas(self, stas):
        self.stas = stas

    def set_freqs(self, freqs):
        self.freqs = freqs

    def set_all(self, M=None, F=None, R=None, N_ev=None, N_sta=None, 
                N_freq=None, data=None, weights=None, scale=None, orids=None,
                stas=None, freqs=None):
        self.M = M
        self.F = F
        self.R = R
        self.N_ev = N_ev
        self.N_sta = N_sta
        self.N_freq = N_freq
        self.data = data
        self.weights = weights
        self.scale = scale
        self.orids = orids
        self.stas = stas
        self.freqs = freqs



class Params(object):
    """
    Parameters class.
    The class features are:
      alpha <--> seismic_moment ln(M_0) [cm^2 * kg^2 / s^2]
      beta <--> corner_frequency  [Hz]
      gamma <--> exponents of geometrical_spreading function 
                                 [adimensional]
      delta <--> quality_factor Q_0 [adimensional] (if Q_alpha is not used)
      site_fi <--> amplitude_correction ln(A) [adimensional]
      site_k  <--> k attenuation factor [s]
      eps_source <--> uncertainty on source contribution
      eps_path <--> uncertainty on path contribution
      eps_site <--> uncertainty on site contribution
      vary <--> flag to fix/release each parameter
      min  <--> lower bound for parameter space
      max  <--> upper bound for parameter space
    """    
    def __init__(self):

        self.alpha = np.array([])
        self.beta = np.array([])
        self.gamma = np.array([])
        self.delta = np.array([])
        self.site_fi = np.array([])
        self.site_k = np.array([])
        self.eps_source = np.array([])
        self.eps_path = np.array([])
        self.eps_site = np.array([])
        self.update_p()
        # self.vary_all(True)

    def update_p(self):
        self._pars = np.concatenate((self.alpha, self.beta,
                                     self.gamma,
                                     self.delta, 
                                     self.site_fi, self.site_k,
                                     self.eps_source, self.eps_path,
                                     self.eps_site))

    @property
    def pars(self):
        self.update_p()
        return self._pars

    def set_pars(self, alpha, beta, gamma, delta, site_fi, site_k,
                 eps_source, eps_path, eps_site):
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.site_fi = site_fi
        self.site_k = site_k
        self.eps_source = eps_source
        self.eps_path = eps_path
        self.eps_site = eps_site
        self.update_p()

    def vary_all(self, value):
        self.vary = np.full(len(self.pars), value)
        if value is True:
            self.set_min_all(-1*np.inf)
            self.set_max_all(1*np.inf)

        if value is False:
            for i in range(len(self.pars)):
                self.min[i] = self.pars[i]-1.e-10
                self.max[i] = self.pars[i]+1.e-10

    # NB 'index' must be a list of indexes to change
    def set_vary(self, index, value):
        for i in range(len(index)):
            self.vary[index[i]] = value
            if value is True:
                self.set_min((index[i],), -1*np.inf)
                self.set_max((index[i],), 1*np.inf)

            if value is False:
                self.set_min((index[i],), (self.pars[index[i]]-1.e-10))
                self.set_max((index[i],), (self.pars[index[i]]+1.e-10))

    def set_min(self, index, value):
        if isinstance(value, np.float64) or np.isinf(value):
            for i in range(len(index)):
                self.min[index[i]] = value
        else:
            for i in range(len(index)):
                self.min[index[i]] = value[i]

    def set_max(self, index, value):
        if isinstance(value, np.float64) or np.isinf(value):
            for i in range(len(index)):
                self.max[index[i]] = value
        else:
            for i in range(len(index)):
                self.max[index[i]] = value[i]

    def set_min_all(self, value):
        self.min = np.full(len(self.pars), value)

    def set_max_all(self, value):
        self.max = np.full(len(self.pars), value)
