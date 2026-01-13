
# =============================================================================
#
# Place for Copyright and authors information
#
# =============================================================================
"""
Implementation of the forward modelling for Fourier 
amplitude spectra for seismic traces, based on the formulas by
Bora et al. (2017), Franceschina et al. (2006) and Edwards et al. (2008)

NB: frequency-dependent site contributions a(f) are not modeled and
left to be calculated from the residuals
"""

import numpy as np
from numba import njit
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils

MOTION = "velocity"
MOTION_TYPE_MAP = {'displacement': 0, 'velocity': 1, 'acceleration': 2}
MOTION_TYPE_INT = MOTION_TYPE_MAP.get(MOTION, 1)

@njit(cache=True)
def scale_to_motion(motion_type, F):
    
    if motion_type == 1:  
        return np.log(2 * _c.PI * F)
    elif motion_type == 2:  
        return 2 * np.log(2 * _c.PI * F)
    return np.zeros_like(F)  

@njit(cache=True)
def ln_q_attenuation_numba(delta_in, F, R):
    
    Q_0, Q_eta = delta_in[0], delta_in[1]
    return -(_c.QCOST / Q_0) * (F**(1. - Q_eta)) * R

def ln_fas(p_in, fcn_args, diff=False):
 
    stat_p, _, p_fixed, index_f, N_gamma, N_delta = fcn_args
    utils.checkext_statparams(stat_p)

    F, R = stat_p.F, stat_p.R
    N_ev, N_sta, N_freq = stat_p.N_ev, stat_p.N_sta, stat_p.N_freq
    scale = stat_p.scale

    p_values = p_in if isinstance(p_in, np.ndarray) else p_in.pars
    p_input = utils.insert_elem(p_values, index_f, p_fixed)
    p = utils.array_to_Param(p_input, stat_p, N_gamma, N_delta)
    
    
    alpha = p.alpha[:, np.newaxis, np.newaxis]
    beta = p.beta[:, np.newaxis, np.newaxis]
    site_fi = p.site_fi[np.newaxis, :, np.newaxis]
    site_k = p.site_k[np.newaxis, :, np.newaxis]
    eps_source = p.eps_source[:, np.newaxis, np.newaxis]
    eps_site = p.eps_site[np.newaxis, :, np.newaxis]
    
    
    F_3d = F.reshape(N_ev, N_sta, N_freq)
    R_3d = R.reshape(N_ev, N_sta, N_freq)

    
    source_spec = -np.log(1 + (F_3d / beta)**2)
    source_term = _c.LN_MCOST + alpha + source_spec
    
    path_geom = utils.ln_piecew_func(R, F, p.gamma, der=None).reshape(N_ev, N_sta, N_freq)
    path_atten = ln_q_attenuation_numba(p.delta, F_3d, R_3d)
    path_term = path_geom + path_atten

    kappa = -_c.PI * (F_3d**(1. - p.delta[1])) * site_k
    site_term = site_fi + kappa
    
    error_term = eps_source + p.eps_path + eps_site
    
    fas = (scale_to_motion(MOTION_TYPE_INT, F_3d) + source_term + path_term + site_term + error_term).flatten()

    if scale is not None:
        fas = fas / scale
    
    if not diff:
        return fas
    else:
        lndata = np.log(stat_p.data)
        if scale is not None:
            lndata = lndata / scale
        return fas - lndata

def handle_func(p_in, *fcn_args):
  
    stat_p, index_v, p_fixed, index_f, N_gamma, N_delta = fcn_args
    fcn_args_tuple = (stat_p, index_v, p_fixed, index_f, N_gamma, N_delta)
    
    diff = ln_fas(p_in, fcn_args_tuple, diff=True)
    
    
    diff *= stat_p.M
    
    weights = stat_p.weights if isinstance(stat_p.weights, np.ndarray) else np.ones_like(diff)
    
    cost = (np.square(diff) * weights).sum()
    
    
    num_points_used = np.sum(stat_p.M)
    if num_points_used > 0:
        cost /= num_points_used
        
    return cost
