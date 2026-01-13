import numpy as np
import spectral_modelling.utils.constants as _c
from spectral_modelling.model_and_invert.fas_log import ln_fas
import spectral_modelling.utils.utils as utils
import gc


def handle_func_deriv(p_in, *fcn_args):
    """
    Analytical Jacobian of the handle (cost) function.

    """

    # unpack
    stat_p, index_v, p_fixed, index_f, N_gamma, N_delta = fcn_args
    utils.checkext_statparams(stat_p)

    M = stat_p.M
    F = stat_p.F
    R = stat_p.R
    N_ev = stat_p.N_ev
    N_sta = stat_p.N_sta
    N_freq = stat_p.N_freq

    # param vector
    p_values = p_in.pars if not isinstance(p_in, np.ndarray) else p_in
    p_input = utils.insert_elem(p_values, index_f, p_fixed)

    
    beta = p_input[N_ev : 2 * N_ev]
    beta = np.tile(np.array([beta]).T, (1, N_sta * N_freq)).reshape(-1)

    gamma = p_input[2 * N_ev : 2 * N_ev + N_gamma]
    delta = p_input[2 * N_ev + N_gamma : 2 * N_ev + N_gamma + N_delta] 

    # (site_k)
    kappa = p_input[2 * N_ev + N_gamma + N_delta + N_sta : 2 * N_ev + N_gamma + N_delta + 2 * N_sta]
    kappa = np.tile(np.array([kappa]).T, N_freq).reshape(-1)
    kappa = np.tile(kappa, N_ev)

  
    wresd = np.zeros_like(p_input)
    ee = np.array([])
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)
    diff = M * ln_fas(p_input, fcn_args_ee, diff=True)

    
    weights = stat_p.weights
    if not isinstance(weights, np.ndarray):
        weights = np.ones_like(diff)
    if len(weights) != len(diff):
        raise utils.MinimizeException("Weights array must be of the same lenght as data")

    # alpha
    diff_alpha = diff * weights
    for l in range(N_ev):
        start = N_sta * N_freq * l
        stop  = N_sta * N_freq * (l + 1)
        wresd[l] = 2.0 * diff_alpha[start:stop].sum()

    # beta
    diff_beta = diff * (F**2 / (beta * (beta**2 + F**2))) * weights
    for l in range(N_ev):
        start = N_sta * N_freq * l
        stop  = N_sta * N_freq * (l + 1)
        wresd[N_ev + l] = 4.0 * diff_beta[start:stop].sum()

    # gamma
    diff_gamma = diff * weights
    for l in range(N_gamma):
        diff_gamma_l = diff_gamma * utils.ln_piecew_func(R, F, gamma, der=l)
        wresd[2 * N_ev + l] = 2.0 * diff_gamma_l.sum()

    # delta
    Q0, Qeta = delta[0], delta[1]
    F_pow = F ** (1.0 - Qeta)

    # d cost / d Q0
    
    diff_dQ0 = diff * F_pow * R * weights
    wresd[2 * N_ev + N_gamma] = (2.0 * _c.PI / (_c.V_S * (Q0**2))) * diff_dQ0.sum()

    # d cost / d Qeta
    
    diff_dQeta = diff * F_pow * np.log(F) * (R / (_c.V_S * Q0) + kappa) * weights
    wresd[2 * N_ev + N_gamma + 1] = 2.0 * _c.PI * diff_dQeta.sum()

    # site_fi
    diff_site_fi = diff * weights
    for l in range(N_sta):
        acc = 0.0
        for i in range(N_ev):
            base = N_sta * N_freq * i
           
            acc += diff_site_fi[base + N_freq * l : base + N_freq * (l + 1)].sum()
        wresd[2 * N_ev + N_gamma + N_delta + l] = 2.0 * acc

    # site_k
    diff_site_k = diff * F_pow * weights
    for l in range(N_sta):
        acc = 0.0
        for i in range(N_ev):
            base = N_sta * N_freq * i
            acc += diff_site_k[base + N_freq * l : base + N_freq * (l + 1)].sum()
        wresd[2 * N_ev + N_gamma + N_delta + N_sta + l] = -2.0 * _c.PI * acc

    # epsilons
    diff_eps = diff * weights

    # eps_source
    for l in range(N_ev):
        start = N_sta * N_freq * l
        stop  = N_sta * N_freq * (l + 1)
        wresd[2 * N_ev + N_gamma + N_delta + 2 * N_sta + l] = 2.0 * diff_eps[start:stop].sum()

    # eps_path 
    wresd[2 * N_ev + N_gamma + N_delta + 2 * N_sta + N_ev] = 2.0 * diff_eps.sum()

    # eps_site
    for l in range(N_sta):
        acc = 0.0
        for i in range(N_ev):
            base = N_sta * N_freq * i
            acc += diff_eps[base + N_freq * l : base + N_freq * (l + 1)].sum()
        wresd[2 * N_ev + N_gamma + N_delta + 2 * N_sta + N_ev + 1 + l] = 2.0 * acc


    wresd = wresd / (N_ev * N_sta * N_freq)

 
    wresd = wresd[index_v]

   
    del diff_alpha, diff_beta, diff_gamma, diff_dQ0, diff_dQeta, diff_site_fi, diff_site_k, diff_eps
    gc.collect()

    return wresd
