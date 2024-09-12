import numpy as np
import os
import pickle
from copy import deepcopy
import spectral_modelling.utils.myClass as myC
from spectral_modelling.utils import config as cfg
import spectral_modelling.utils.constants as _c
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing 
import scipy.fftpack as fftpack
from obspy import Stream, Trace
from scipy.interpolate import CubicSpline
from spectral_modelling.db_create.f1f2 import f1f2


# -----------------------------------------------------------------------------
# Read pickle object containing a list
# -----------------------------------------------------------------------------
def read_pickle_list(pklname):
    pkllist = []
    with open(pklname, 'rb') as f:
        while True:
            try:
                pkllist.append(pickle.load(f))
            except EOFError:
                break
    return pkllist


# -----------------------------------------------------------------------------
# Read input component from pickle object
# -----------------------------------------------------------------------------
def read_input_component(runname=cfg.RUNNAME, sptdbpath=cfg.SPTDBPATH,
                         hcomponent='R', component='S'):
    """Read input pickle files:
       hcomponent = ['R', 'T']
         - R = radial horizontal
         - T = transverse horizontal
       component = ['S', 'N']
         - S = signal
         - N = noise

       Returns:
         orids, stas, freqs = lists of origin IDs, station names and 
                              frequency points
         ml         = list of magnitudes (ml) of the events 
         M, F, R, Z = arrays of boolean masks, frequency points,
                      hypocentral distances (in cm) and smoothed
                      spectral amplitudes
         fmin, fmax = arrays of minimum and maximum used frequencies
    """

    if component not in ('S', 'N'):
        raise MinimizeException("Component must be either 'S' = signal or 'N' = noise")
    
    if hcomponent not in ('R', 'T'):
        raise MinimizeException("Horizontal component must be either 'R' = radial or 'T' = transverse")

    pklname = os.path.join(sptdbpath, runname + cfg.COMPNAME[component] + hcomponent + '.pkl')
    w_list = read_pickle_list(pklname)

    orids = []
    stas = []
    ml = []
    mag = []
    for w in w_list:
        if w.orid not in orids:
            orids.append(w.orid)
            ml.append(w.ml)
            mag.append(w.mag)
        if w.sta not in stas:
            stas.append(w.sta)

    ml = np.array(ml)
    mag = np.array(mag)
    freqs = w_list[0].sm_freq

    N_i = len(orids)
    N_j = len(stas)
    N_k = len(freqs)
    orids_dict = dict(zip(orids, list(range(N_i))))
    stas_dict = dict(zip(stas, list(range(N_j))))

    M = np.zeros((N_i * N_j * N_k))
    F = np.full((N_i * N_j * N_k), 1.e-20)
    R = np.full((N_i * N_j * N_k), 1.e-20)
    Z = np.full((N_i * N_j * N_k), 1.e-20)
    fmin = np.zeros(N_i * N_j)
    fmax = np.zeros(N_i * N_j)

    for w in w_list:
        i = orids_dict[w.orid]
        j = stas_dict[w.sta]
        fmin[N_j * i + j] = w.freq_lims[0]
        fmax[N_j * i + j] = w.freq_lims[1]
        for k in range(N_k):
            M[N_j * N_k * i + N_k * j + k] = w.fbool[k]
            F[N_j * N_k * i + N_k * j + k] = w.sm_freq[k]
            R[N_j * N_k * i + N_k * j + k] = w.hypdist * 100000  #to CGS
            Z[N_j * N_k * i + N_k * j + k] = w.sm_amp[k]

    stat_p = myC.StaticParams(M=M, F=F, R=R, N_ev=N_i, N_sta=N_j, N_freq=N_k,
                          data=Z, orids=orids, stas=stas, freqs=freqs)
    return stat_p, ml, mag, fmin, fmax



def getSNR(S, N):
    """
    Function to calculate SNR of 'signal' waveform (S) relative to
    pre-signal noise (N). 
    """
    checkext_trace(S)
    checkext_trace(N)

    if (np.size(S.data) < 1) or (np.size(N.data) < 1):
        print("getSNR: Empty signal or noise data")
        return None

    # demean over entire waveform
    S.data = S.data - np.mean(N.data)
    N.data = N.data - np.mean(N.data)

    # calculate ratios
    signallevel = np.abs(S.data).max()
    noiselevel = np.sqrt(np.mean(np.square(N.data)))

    SNR = signallevel / noiselevel
    SNRdB = 10 * np.log10(SNR)
    return SNRdB


def fas_calc_and_store(stream_s, stream_n, spectrum_s, spectrum_n, index):
    """
    Utility to calculate FAS of given Stream[index] and store it in 
    the given Spectrum object, for corresponding signal and noise streams
    """
    checkext_stream(stream_s)
    checkext_stream(stream_n)
    checkext_spectrum(spectrum_s)
    checkext_spectrum(spectrum_n)
    freqs = cfg.FREQS

    # 1 - signal
    tr_s = stream_s[index]
    npts_s = tr_s.stats["npts"]
    dt_s = tr_s.stats["delta"]
    # Calculate amplitude fft spectra
    nfrqs_s = int(npts_s/2)+1 if (npts_s % 2) == 0 else int((npts_s+1)/2)
    frq_s = np.abs(np.fft.fftfreq(npts_s, dt_s))[:nfrqs_s]
    hw_s = np.array([(1 - np.cos(2 * _c.PI * j / npts_s))/2. 
                    for j in range(npts_s)])
    tr_s_fft = np.sqrt(2*((np.abs(fftpack.fft(tr_s.data)[:nfrqs_s]))**2)
                        / (np.sum(hw_s**2)) / dt_s)
    # Calculate KO smoothed fft and resample it over FREQS
    sm_s_fft = konno_ohmachi_smoothing(tr_s_fft, frq_s, bandwidth=cfg.KO_B,
                                     normalize=True)
    cs_s_sm = CubicSpline(frq_s, sm_s_fft)

    # 2 - noise
    tr_n = stream_n[index]
    npts_n = tr_n.stats["npts"]
    dt_n = tr_n.stats["delta"]
    # Calculate amplitude fft spectra
    nfrqs_n = int(npts_n/2)+1 if (npts_n % 2) == 0 else int((npts_n+1)/2)
    frq_n = np.abs(np.fft.fftfreq(npts_n, dt_n))[:nfrqs_n]
    hw_n = np.array([(1 - np.cos(2 * _c.PI * j / npts_n))/2. 
                    for j in range(npts_n)])
    tr_n_fft = np.sqrt(2*((np.abs(fftpack.fft(tr_n.data)[:nfrqs_n]))**2)
                        / (np.sum(hw_n**2)) / dt_n)
    # Calculate KO smoothed fft and resample it over FREQS
    sm_n_fft = konno_ohmachi_smoothing(tr_n_fft, frq_n, bandwidth=cfg.KO_B,
                                     normalize=True)
    cs_n_sm = CubicSpline(frq_n, sm_n_fft)

    flims = f1f2(cs_s_sm(freqs), cs_n_sm(freqs), freqs)
    if flims is not None:
        fbools = np.zeros_like(freqs)
        index1 = [i for i in range(len(freqs)) if (freqs[i] >= flims[0]) and (freqs[i] <= flims[1])]
        fbools[index1] = 1

        spectrum_s.freq = frq_s
        spectrum_s.amp = tr_s_fft
        spectrum_s.sm_freq = freqs
        spectrum_s.sm_amp = cs_s_sm(freqs)
        spectrum_s.fbool = fbools 
        spectrum_s.freq_lims = flims

        spectrum_n.freq = frq_n
        spectrum_n.amp = tr_n_fft
        spectrum_n.sm_freq = freqs
        spectrum_n.sm_amp = cs_n_sm(freqs)
        spectrum_n.fbool = fbools 
        spectrum_n.freq_lims = flims
        return 1

    return 0

# -----------------------------------------------------------------------------
# Utility for re-building whole Params object  (from fixed and free
# arrays)
# -----------------------------------------------------------------------------
def insert_elem(a, insert_at, values):
    b = list(a[:])
    insert_at = list(insert_at)
    values = list(values)

    for i in range(len(insert_at)):
        b[insert_at[i]:insert_at[i]] = [values[i]]

    return np.array(b)



# -----------------------------------------------------------------------------
# Piecewise function utility used for geometrical spreading term - Malagnini
# -----------------------------------------------------------------------------
def ln_piecew_func(x, f, g, der=None):
    """Geometrical spreading piecewise function and derivative (up to 5
       intervals).
       x = hypocentral distances
       f = corresponding used frequencies
       g = exponents to be used
       der = derivative status (None=no, i = index of g component
                                respect to which we derive)
    """

    index_lf = np.where(f <= 1.)[0]
    index_hf = np.where(f > 1.)[0]

    x_lf = x[index_lf]
    conditions_lf = [(x_lf <= _c.R2),
                     ((x_lf > _c.R2) & (x_lf <= _c.R3)),
                     ((x_lf > _c.R3) & (x_lf <= _c.R4)),
                     ((x_lf > _c.R4) & (x_lf <= _c.R5)),
                     (x_lf > _c.R5)]
    choices_lf = [g[0]*np.log(x_lf/_c.R0),
                  g[1]*np.log(x_lf/_c.R2) + g[0]*np.log(_c.R2/_c.R0),
                  g[2]*np.log(x_lf/_c.R3) + g[1]*np.log(_c.R3/_c.R2)
                  + g[0]*np.log(_c.R2/_c.R0),
                  g[3]*np.log(x_lf/_c.R4) + g[2]*np.log(_c.R4/_c.R3)
                  + g[1]*np.log(_c.R3/_c.R2)
                  + g[0] * np.log(_c.R2/_c.R0),
                  g[4]*np.log(x_lf/_c.R5) + g[3]*np.log(_c.R5/_c.R4)
                  + g[2]*np.log(_c.R4/_c.R3) + g[1]*np.log(_c.R3/_c.R2)
                  + g[0]*np.log(_c.R2/_c.R0)]

    x_hf = x[index_hf]
    conditions_hf = [(x_hf <= _c.R1),
                     ((x_hf > _c.R1) & (x_hf <= _c.R2)),
                     ((x_hf > _c.R2) & (x_hf <= _c.R3)),
                     ((x_hf > _c.R3) & (x_hf <= _c.R5)),
                     (x_hf > _c.R5)]
    choices_hf = [g[5]*np.log(x_hf/_c.R0),
                  g[6]*np.log(x_hf/_c.R1) + g[5]*np.log(_c.R1/_c.R0),
                  g[7]*np.log(x_hf/_c.R2) + g[6]*np.log(_c.R2/_c.R1)
                  + g[5]*np.log(_c.R1/_c.R0),
                  g[8]*np.log(x_hf/_c.R3) + g[7]*np.log(_c.R3/_c.R2)
                  + g[6]*np.log(_c.R2/_c.R1)
                  + g[5]*np.log(_c.R1/_c.R0),
                  g[9]*np.log(x_hf/_c.R5) + g[8]*np.log(_c.R5/_c.R3)
                  + g[7]*np.log(_c.R3/_c.R2) + g[6]*np.log(_c.R2/_c.R1)
                  + g[5]*np.log(_c.R1/_c.R0)]

    if der == 0:
        choices_lf = [np.log(x_lf/_c.R0), np.log(_c.R2/_c.R0),
                      np.log(_c.R2/_c.R0), np.log(_c.R2/_c.R0),
                      np.log(_c.R2/_c.R0)]
        choices_hf = [0, 0, 0, 0, 0]
    if der == 1:
        choices_lf = [0, np.log(x_lf/_c.R2), np.log(_c.R3/_c.R2),
                      np.log(_c.R3/_c.R2), np.log(_c.R3/_c.R2)]
        choices_hf = [0, 0, 0, 0, 0]
    if der == 2:
        choices_lf = [0, 0, np.log(x_lf/_c.R3), np.log(_c.R4/_c.R3),
                      np.log(_c.R4/_c.R3)]
        choices_hf = [0, 0, 0, 0, 0]
    if der == 3:
        choices_lf = [0, 0, 0, np.log(x_lf/_c.R4), np.log(_c.R5/_c.R4)]
        choices_hf = [0, 0, 0, 0, 0]
    if der == 4:
        choices_lf = [0, 0, 0, 0, np.log(x_lf/_c.R5)]
        choices_hf = [0, 0, 0, 0, 0]
    if der == 5:
        choices_lf = [0, 0, 0, 0, 0]
        choices_hf = [np.log(x_hf/_c.R0), np.log(_c.R1/_c.R0),
                      np.log(_c.R1/_c.R0), np.log(_c.R1/_c.R0),
                      np.log(_c.R1/_c.R0)]
    if der == 6:
        choices_lf = [0, 0, 0, 0, 0]
        choices_hf = [0, np.log(x_hf/_c.R1), np.log(_c.R2/_c.R1),
                      np.log(_c.R2/_c.R1), np.log(_c.R2/_c.R1)]
    if der == 7:
        choices_lf = [0, 0, 0, 0, 0]
        choices_hf = [0, 0, np.log(x_hf/_c.R2), np.log(_c.R3/_c.R2),
                      np.log(_c.R3/_c.R2)]
    if der == 8:
        choices_lf = [0, 0, 0, 0, 0]
        choices_hf = [0, 0, 0, np.log(x_hf/_c.R3), np.log(_c.R5/_c.R3)]
    if der == 9:
        choices_lf = [0, 0, 0, 0, 0]
        choices_hf = [0, 0, 0, 0, np.log(x_hf/_c.R5)]

    pieces_lf = np.select(conditions_lf, choices_lf)
    pieces_hf = np.select(conditions_hf, choices_hf)
    pieces = np.zeros_like(x)
    for i in range(len(pieces_lf)):
        pieces[index_lf[i]] = pieces_lf[i]
    for i in range(len(pieces_hf)):
        pieces[index_hf[i]] = pieces_hf[i]
    return pieces



# -----------------------------------------------------------------------------
# Select stations to be used as reference
# -----------------------------------------------------------------------------
def stas_flat(runname=cfg.RUNNAME, sptdbpath=cfg.SPTDBPATH,
              hcomponent='R'):
    """Create condition on reference stations.
       Defaults to constraining the average amplification of sites on
       EC8 class A to zero. The stations to be included in the
       constraint are read from the soil_dict dictionary inside the
       constants utility.
       Returns array with indices of stations to be used in the reference
    """

    stat_p, ml, fmin, fmax, hcompout = \
        read_input_component(runname=runname, sptdbpath=sptdbpath,
                             hcomponent=hcomponent, component='S')
    stas = stat_p.stas
    N_j = stat_p.N_sta
    stas_dict = dict(zip(stas, list(range(N_j))))

    # meanA
    flatstas = []
    for s in range(len(_c.soil_dict)):
        if list(_c.soil_dict.values())[s] == 'A':
            flatstas.append(list(_c.soil_dict)[s])

    flatindex = np.array([stas_dict.get(x) for x in flatstas])

    return flatindex



# -----------------------------------------------------------------------------
# Build input pars with free terms only and additional array of fixed
# terms
# -----------------------------------------------------------------------------
def create_p_in(p_input, stat_p):
    """Create p_in from p_input, keeping only pars free to vary.
       NOTE p_input has to be a Params() instance
       Returns:
         p_in    = Params array of free parameter values
         index_v = array of indexes marking the original position of 
                   free parameters
         p_fixed = array of fixed paramater values
         index_f = array of indexes marking the original position of 
                   fixed parameters
    """

    checkext(p_input, stat_p)
    N_ev = stat_p.N_ev
    N_sta = stat_p.N_sta
    N_gamma = p_input.gamma.size
    N_delta = p_input.delta.size

    # read input parameters and associated status (fixed or not)
    mask = p_input.vary
    if len(mask) == 0:
        mask = np.ones_like(p_input.pars)

    # define indices for slicing
    index_v = np.where(mask == 1)[0]
    index_f = np.where(mask == 0)[0]

    idx_v_a = index_v[np.where(index_v < N_ev)]
    idx_v_b = index_v[np.where((N_ev <= index_v) & (index_v < 2*N_ev))]
    idx_v_c = index_v[np.where((2*N_ev <= index_v) &
                                (index_v < (2*N_ev + N_gamma)))]
    idx_v_d = index_v[np.where(((2*N_ev + N_gamma) <= index_v) &
                                (index_v < (2*N_ev + N_gamma + N_delta)))]
    idx_v_sfi = index_v[np.where(((2*N_ev + N_gamma + N_delta) <= index_v) &
                                  (index_v < (2*N_ev + N_gamma + N_delta + N_sta)))]
    idx_v_sk = index_v[np.where(((2*N_ev + N_gamma + N_delta + N_sta) <= index_v) &
                                 (index_v < (2*N_ev + N_gamma + N_delta + 2*N_sta)))]
    idx_v_eso = index_v[np.where(((2*N_ev + N_gamma + N_delta + 2*N_sta) <= index_v) &
                                  (index_v < (2*N_ev + N_gamma + N_delta + 2*N_sta + N_ev)))]
    idx_v_epa = index_v[np.where(((2*N_ev + N_gamma + N_delta + 2*N_sta + N_ev) <= index_v) &
                                  (index_v < (2*N_ev + N_gamma + N_delta + 2*N_sta + N_ev + 1)))]
    idx_v_esi = index_v[np.where((2*N_ev + N_gamma + N_delta + 2*N_sta + N_ev + 1) <= index_v)]

    idx_v_b = idx_v_b - N_ev
    idx_v_c = idx_v_c - 2*N_ev
    idx_v_d = idx_v_d - (2*N_ev + N_gamma)
    idx_v_sfi = idx_v_sfi - (2*N_ev + N_gamma + N_delta)
    idx_v_sk = idx_v_sk - (2*N_ev + N_gamma + N_delta + N_sta)
    idx_v_eso = idx_v_eso - (2*N_ev + N_gamma + N_delta + 2*N_sta)
    idx_v_epa = idx_v_epa - (2*N_ev + N_gamma + N_delta + 2*N_sta + N_ev)
    idx_v_esi = idx_v_esi - (2*N_ev + N_gamma + N_delta + 2*N_sta + N_ev + 1)

    # slice
    alpha_in = p_input.alpha[idx_v_a]
    beta_in = p_input.beta[idx_v_b]
    gamma_in = p_input.gamma[idx_v_c]
    delta_in = p_input.delta[idx_v_d]
    site_fi_in = p_input.site_fi[idx_v_sfi]
    site_k_in = p_input.site_k[idx_v_sk]
    eps_source_in = p_input.eps_source[idx_v_eso]
    eps_path_in = p_input.eps_path[idx_v_epa]
    eps_site_in = p_input.eps_site[idx_v_esi]

    # build new parameter object with variable parameters

    p_in = myC.Params()
    p_in.set_pars(alpha_in, beta_in, gamma_in, delta_in, site_fi_in,
                  site_k_in, eps_source_in, eps_path_in, eps_site_in)

    # build new array object with fixed parameters
    p_fixed = p_input.pars[index_f]

    return p_in, index_v, p_fixed, index_f


# -----------------------------------------------------------------------------
# Utility for creating a stripped version of a Params object  
# (same parameters except for the three epsilon terms, which are set to
#  zero)
# (stat_p is used to retrieve correct dimensions only)
# NB!! Input and output are Param objects
# -----------------------------------------------------------------------------
def strip_pars(p_in, stat_p):
    checkext(p_in, stat_p)
    p_strip = deepcopy(p_in)

    eps_source_strip = np.zeros(stat_p.N_ev)
    eps_path_strip = np.array([0.])
    eps_site_strip = np.zeros(stat_p.N_sta)

    p_strip.set_pars(p_in.alpha, p_in.beta, p_in.gamma, p_in.delta,
                     p_in.site_fi, p_in.site_k, eps_source_strip,
                     eps_path_strip, eps_site_strip)

    return p_strip


# -----------------------------------------------------------------------------
# Utility for transforming a parameter array into a Params object  
# (stat_p is used to retrieve correct dimensions only)
# -----------------------------------------------------------------------------
def array_to_Param(p_in, stat_p, N_gamma, N_delta):
    # safecheck in case p_in is not an array
    checkext_ndarray(p_in)

    N_i = stat_p.N_ev
    N_j = stat_p.N_sta

    alpha = p_in[0:N_i]
    beta = p_in[N_i:(2*N_i)]
    gamma = p_in[(2*N_i):(2*N_i + N_gamma)]
    delta = p_in[(2*N_i + N_gamma):(2*N_i + N_gamma + N_delta)]
    site_fi = p_in[(2*N_i + N_gamma + N_delta):
                   (2*N_i + N_gamma + N_delta + N_j)]
    site_k = p_in[(2*N_i + N_gamma + N_delta + N_j):
                  (2*N_i + N_gamma + N_delta + 2*N_j)]
    eps_source = p_in[(2*N_i + N_gamma + N_delta + 2*N_j):
                      (2*N_i + N_gamma + N_delta + 2*N_j + N_i)]
    eps_path = p_in[(2*N_i + N_gamma + N_delta + 2*N_j + N_i):
                    (2*N_i + N_gamma + N_delta + 2*N_j + N_i + 1)]
    eps_site = p_in[(2*N_i + N_gamma + N_delta + 2*N_j + N_i + 1):
                    (2*N_i + N_gamma + N_delta + 2*N_j + N_i + 1 + N_j)]

    p_out = myC.Params()
    p_out.set_pars(alpha, beta, gamma, delta, site_fi, site_k, eps_source,
                  eps_path, eps_site)

    return p_out



def fcmin_NE_IT(m):
    """
    Utility to calculate lower bound on corner frequency for
    North-East Italy, using literature regional spans for stress drop
    and seismic moment
    """
    fc_min = m0_to_beta_Brune(m0=5.623 * m, stdrop=1.e+06)
    return fc_min


def fcmax_NE_IT(m):
    """
    Utility to calculate upper bound on corner frequency for
    North-East Italy, using literature regional spans for stress drop
    and seismic moment
    """
    fc_max = m0_to_beta_Brune(m0=0.178 * m, stdrop=5.e+07)
    return fc_max



# -----------------------------------------------------------------------------
# Exception handler utilities
# -----------------------------------------------------------------------------
class MinimizeException(Exception):
    """General Purpose Exception."""

    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg

    def __str__(self):
        return "{}".format(self.msg)


def checkext_spectrum(spt):
    """Check extension to see if a Spectrum object."""
    if not isinstance(spt, myC.Spectrum):
        raise MinimizeException(
            "spectrum must be a myClass.Spectrum() instance")

def checkext_stream(st):
    """Check extension to see if a obspy Stream object."""
    if not isinstance(st, Stream):
        raise MinimizeException(
            "stream must be an obspy Stream() instance")

def checkext_trace(tr):
    """Check extension to see if a obspy Stream object."""
    if not isinstance(tr, Trace):
        raise MinimizeException(
            "trace must be an obspy Trace() instance")

def checkext(params, stat_params):
    """Check extension to see if a Params/StaticParams object."""

    if not isinstance(params, myC.Params):
        raise MinimizeException(
            "params (p_in) must be a myClass.Params() instance")

    if not isinstance(stat_params, myC.StaticParams):
        raise MinimizeException(
            "stat_params (stat_p) must be a myClass.StaticParams() instance")


def checkext_params(params):
    """Check extension to see if a Params object."""
    if not isinstance(params, myC.Params):
        raise MinimizeException(
            "params (p_in) must be a myClass.Params() instance")


def checkext_statparams(stat_params):
    """Check extension to see if a StaticParams object."""
    if not isinstance(stat_params, myC.StaticParams):
        raise MinimizeException(
            "stat_params (stat_p) must be a myClass.StaticParams() instance")


def checkext_ndarray(arr):
    """Check extension to see if a np.ndarray object."""
    if not isinstance(arr, np.ndarray):
        raise MinimizeException("p_arr must be a np.ndarray instance")

# -----------------------------------------------------------------------------
# Converter utilities used in all calculations
# -----------------------------------------------------------------------------
def ml_to_mag_munafo(ml):
    """
    Converter from local magnitude (ml) to moment magnitude (mag) 
    NOTE: it uses formula from Munafo' et al. (2016) for small Italian
    earthquakes 
    """
    mag = 2. * ml / 3. + 1.15
    return mag


def mag_to_m0_HK(mag):
    """
    Converter from moment magnitude (mag) to seismic moment 
    (m0 [dyn* cm = g*cm^2/s^2]), using formula from Hanks and Kanamori
    (1979)
    """
    m0 = 10 ** (1.5 * mag + 16.1)
    return m0


def m0_to_beta_Brune(m0, stdrop):
    """
    Converter from seismic moment (m0) to beta (corner frequency [Hz]), 
    using formula from Brune (1970, 1971) source model
    """
    beta = 0.4906 * _c.V_S * np.power(stdrop / m0, (1. / 3.))
    return beta

def mag_to_ml_munafo(mag):
    """
    Converter from moment magnitude (mag) to local magnitude (ml) 
    NOTE: it inverts formula from Munafo' et al. (2016) for small Italian
    earthquakes 
    """
    ml = (mag - 1.15) * 3. / 2.
    return ml

def stress_drop_Brune(m0, beta):
    """
    Calculate stress drop from seismic moment (m0) and beta 
    (corner frequency [Hz]), using formula from Brune (1970, 1971) source model
    """
    stdrop = ((beta/(0.4906*_c.V_S))**3) * m0
    return stdrop

def m0_to_mag_HK(m0):
    """
    Converter from seismic moment (m0 [dyn* cm = g*cm^2/s^2]) to moment 
    magnitude (mag), using formula from Hanks and Kanamori (1979)
    """
    mag = 2. * (np.log10(m0) - 16.1) / 3.
    return mag


# -----------------------------------------------------------------------------
# Logging utilities
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Utility for building the inversion run label from choices made on used
# model, uncertainty handling, weights, reference amplitude, bounds,
# noise, waveform database version and inversion method
# -----------------------------------------------------------------------------
def set_run_label(model_name, use_uncert, weights, bounds, subtract_noise,
                  ref_ampl, method, hcomponent):
    noisestatus = 'nonoise_'
    if subtract_noise == True:
        noisestatus = 'noise_'
    weightstatus = '_noweights_'
    if weights is not None:
        weightstatus = '_weighted_'
    boundstatus = '_nobounds_'
    if bounds is not None:
        boundstatus = '_bounded_'

    label = use_uncert + '_' + model_name + weightstatus + ref_ampl \
            + boundstatus + noisestatus + method + '_' + hcomponent
    return label


# -----------------------------------------------------------------------------
# Utility for reading labels regarding inversion choices from the
# inversion run label
# -----------------------------------------------------------------------------
def read_run_label(label):
    use_uncert = label.split('_')[0]
    model_name = label.split('_')[1]
    weightstatus = label.split('_')[2]
    ref_ampl = label.split('_')[3]
    boundstatus = label.split('_')[4]
    noisestatus = label.split('_')[5]
    method = label.split('_')[6]
    hcomponent = label.split('_')[7]

    subtract_noise = 'False'
    if noisestatus == 'noise':
        subtract_noise = 'True'
    weights = 'False'
    if weightstatus == 'weighted':
        weights = 'True'
    bounds = 'False'
    if boundstatus == 'bounded':
        bounds = 'True'

    choices = use_uncert, model_name, weights, ref_ampl, bounds, \
              subtract_noise, method, hcomponent
    return choices


def create_pars_dict(params, stat_p):
    """Create dictionary of parameters"""

    checkext_params(params)
    index = []
    names = []
    N_i = stat_p.N_ev
    N_j = stat_p.N_sta
    N_gamma = params.gamma.size
    N_delta = params.delta.size
    for i in range(len(params.pars)):
        index.append(i)
        if i < N_i:
            names.append('alpha_' + str(i))
        if N_i <= i < 2*N_i:
            names.append('beta_' + str(i - N_i))
        if 2*N_i <= i < (2*N_i + N_gamma):
            names.append('gamma_' + str(i - (2*N_i)))
        if (2*N_i + N_gamma) <= i < (2*N_i + N_gamma + N_delta):
            names.append('delta_' + str(i - (2*N_i + N_gamma)))
        if (2*N_i + N_gamma + N_delta) <= i \
                < (2*N_i + N_gamma + N_delta + N_j):
            names.append('site_fi_' + str(i - (2*N_i + N_gamma + N_delta)))
        if (2*N_i + N_gamma + N_delta + N_j) <= i \
                < (2*N_i + N_gamma + N_delta + 2*N_j):
            names.append('site_k_' + str(i - (2*N_i + N_gamma + N_delta
                                              + N_j)))
        if (2*N_i + N_gamma + N_delta + 2*N_j) <= i \
                < (2*N_i + N_gamma + N_delta + 2*N_j + N_i):
            names.append('eps_source_' + str(i - (2*N_i + N_gamma + N_delta
                                                  + 2*N_j)))
        if (2*N_i + N_gamma + N_delta + 2*N_j + N_i) <= i \
                < (2*N_i + N_gamma + N_delta + 2*N_j + N_i + 1):
            names.append('eps_path_' + str(i - (2*N_i + N_gamma + N_delta
                                                + 2*N_j + N_i)))
        if i >= (2*N_i + N_gamma + N_delta + 2*N_j + N_i + 1):
            names.append('eps_site_' + str(i - (2*N_i + N_gamma + N_delta
                                                + 2*N_j + N_i + 1)))

    names = np.array(names)
    index = np.array(index)
    pars_dict = dict(zip(index, names))
    return pars_dict


def writelog(logpath, filename, output, pars_dict, p_in, p_true, params):
    """Writes log of output parameters of the inversion.
       Note that p_true is only useful for synthetic tests, 
       otherwise it can be set equal to p_in
    """
    if os.path.exists(logpath) == 0:
        os.makedirs(logpath)
    logfile = open(logpath + '/' + filename, 'w')
    logfile.write('#status: ' + str(output.status) + '\n')
    logfile.write('#success: ' + str(output.success) + '\n')
    logfile.write('#message: ' + str(output.message) + '\n')
    logfile.write('#fun : ' + str(output.fun) + '\n')
    try:
        logfile.write('#nfev, njev: ' + str(output.nfev) + ' '
                      + str(output.njev) + '\n' + '#\n')
    finally:
        logfile.write('#nfev: ' + str(output.nfev) + '\n' + '#\n')
    logfile.write('{:12}'.format('#parameter') + '{:24}'.format('input value')
                  + '{:24}'.format('true value')
                  + '{:26}'.format('output value')
                  + '{:10}'.format('diff. (%)') + '\n')
    for i in range(len(params.pars)):
        if p_true.pars[i] != 0:
            percent_diff = (100 * (p_true.pars[i] - params.pars[i])
                            / p_true.pars[i])
            logfile.write('{:18}'.format(str(pars_dict[i]))
                          + '{:24}'.format(str(p_in.pars[i]))
                          + '{:24}'.format(str(p_true.pars[i]))
                          + '{:26}'.format(str(params.pars[i]))
                          + '{:10}'.format(str('%.4f' % percent_diff)) + '\n')
        else:
            logfile.write('{:18}'.format(str(pars_dict[i]))
                          + '{:24}'.format(str(p_in.pars[i]))
                          + '{:24}'.format(str(p_true.pars[i]))
                          + '{:26}'.format(str(params.pars[i])) + 'n.d.'
                          + '\n')
    logfile.close()

