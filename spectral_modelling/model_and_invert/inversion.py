import numpy as np
import os
import pickle
import timeit
from copy import deepcopy
from tqdm import tqdm

import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils
import spectral_modelling.utils.myClass as myC
from spectral_modelling.model_and_invert.fas_log import handle_func
from spectral_modelling.model_and_invert.jacobian_module import handle_func_deriv
from spectral_modelling.model_and_invert.wrap_function import minimizer_wrapper


def bounds_NE_Italy(pars, alpha_ref, m0_ref):
    """  Define initial parameter bounds for seismic setting of
    North East Italy"""
    #alpha
    lb_alpha = alpha_ref - 1.727
    ub_alpha = alpha_ref + 1.727
    
    #beta
    lb_beta = utils.fcmin_NE_IT(m0_ref)
    ub_beta = utils.fcmax_NE_IT(m0_ref)
    
    #gamma
    lb_gamma = np.full_like(pars.gamma, -2.0)
    ub_gamma = np.full_like(pars.gamma, -0.3)
    
    #delta
    lb_delta = np.array([10., 0.])
    ub_delta = np.array([500., 1.])
    
    #site_fi and site_k
    lb_site_fi = np.full_like(pars.site_fi, -1.61)
    ub_site_fi = np.full_like(pars.site_fi, 1.61)
    lb_site_k = np.full_like(pars.site_k, 0.01)
    ub_site_k = np.full_like(pars.site_k, 0.15)
    
    #epsilons
    num_eps = len(pars.eps_source) + len(pars.eps_path) + len(pars.eps_site)
    lb_eps = np.full(num_eps, -np.inf)
    ub_eps = np.full(num_eps, np.inf)

    
    lb = np.concatenate([lb_alpha, lb_beta, lb_gamma, lb_delta, lb_site_fi, lb_site_k, lb_eps])
    ub = np.concatenate([ub_alpha, ub_beta, ub_gamma, ub_delta, ub_site_fi, ub_site_k, ub_eps])
    
    return (lb, ub)


def fix_params(pars, model_name, use_uncert, ref_ampl, stas_dict):
    if len(pars.vary) == 0:
        pars.vary_all(True)
    N_i, N_j = pars.alpha.size, pars.site_fi.size
    N_gamma, N_delta = pars.gamma.size, pars.delta.size
    if model_name == 'malagniniQ0':
        pars.set_vary((2*N_i + N_gamma + 1,), False)
    if model_name == 'gamma-deltas':
        for i in range(N_delta):
            pars.set_vary((2*N_i + N_gamma + i,), False)
    if model_name.startswith('fixA'):
        for j in range(N_j):
            pars.set_vary((2*N_i + N_gamma + N_delta + j,), False)
    if use_uncert == '2eps':
        for j in range(N_j):
            pars.set_vary((2*N_i + N_gamma + N_delta + 2*N_j + N_i + 1 + j,), False)
    if use_uncert == 'noeps':
        for j in range(N_j + N_i + 1):
            pars.set_vary((2*N_i + N_gamma + N_delta + 2*N_j + j,), False)
    if ref_ampl == 'PURA':
        pars.set_vary((2*N_i + N_gamma + N_delta + stas_dict['PURA'],), False)
    return pars


class HandleFuncWithProgress:
    def __init__(self, fcn, stat_p, index_v, p_fixed, index_f, N_gamma, N_delta, total_calls=10000):
        self.fcn = fcn
        self.fcn_args = (stat_p, index_v, p_fixed, index_f, N_gamma, N_delta)
        self.pbar = tqdm(total=total_calls, desc="Inversion progress", ncols=90)
        self.calls = 0

    def __call__(self, p_in, *args):
        self.calls += 1
        self.pbar.update(1)
        return self.fcn(p_in, *args)

    def close(self):
        self.pbar.close()



def fas_invert(runname=cfg.RUNNAME, sptdbpath=cfg.SPTDBPATH, outpath=cfg.RUNPATH,
               hcomponent='R', model_name='malagniniQ0', method='SLSQP',
               jac='exact', use_uncert='noeps', ref_ampl='meanA',
               subtract_noise=False, weights=None, bounds='NE_Italy', save=False):

    DEFAULT_OPTIONS = {'disp': True, 'maxiter': 10000}

    stat_p_s, ml_s, fmin_s, fmax_s, _ = utils.read_input_component(
        runname=runname, sptdbpath=sptdbpath, hcomponent=hcomponent, component='S')

    N_i, N_j, N_k = stat_p_s.N_ev, stat_p_s.N_sta, stat_p_s.N_freq
    stas_dict = dict(zip(stat_p_s.stas, list(range(N_j))))
    mag = utils.ml_to_mag_munafo(ml_s)
    m0_calc = utils.mag_to_m0_HK(mag)
    alpha_calc = np.log(m0_calc)

    with open(f"{cfg.SYNTHPARS_PATH}/{runname}_{model_name}.pkl", 'rb') as f:
        p_model = pickle.load(f)

    N_gamma, N_delta = p_model.gamma.size, p_model.delta.size
    run_name = utils.set_run_label(model_name, use_uncert, weights, bounds, subtract_noise, ref_ampl, method, hcomponent)
    runpath = os.path.join(cfg.RUNPATH, run_name)
    os.makedirs(runpath, exist_ok=True)

    p_model = fix_params(p_model, model_name, use_uncert, ref_ampl, stas_dict)
    stat_p = myC.StaticParams(M=stat_p_s.M, F=stat_p_s.F, R=stat_p_s.R,
                              N_ev=N_i, N_sta=N_j, N_freq=N_k)
    pars_dict = utils.create_pars_dict(p_model, stat_p)

    if subtract_noise:
        stat_p_n, _, _, _, _ = utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                                          hcomponent=hcomponent, component='N')
        Z_diff = stat_p_s.data - stat_p_n.data
        Z_true = np.where(Z_diff > 0, Z_diff, 1.0)
    else:
        Z_true = stat_p_s.data

    stat_p.set_data(Z_true)
    p_input = deepcopy(p_model)

    if bounds == 'NE_Italy':
        inv_bounds = bounds_NE_Italy(p_model, alpha_calc, m0_calc)
    elif bounds is not None and len(bounds) != len(p_model.pars):
        raise utils.MinimizeException("Bounds must match length of parameter object.")
    else:
        inv_bounds = bounds

   
    stat_p_check, index_v, p_fixed, index_f = utils.create_p_in(p_input, stat_p)
    fcn_with_progress = HandleFuncWithProgress(handle_func, stat_p, index_v, p_fixed, index_f, N_gamma, N_delta)

    kwargs = dict(method='SLSQP', bounds=inv_bounds, options=DEFAULT_OPTIONS)
    if jac == 'exact':
        kwargs['jac'] = handle_func_deriv
    else:
        kwargs['jac'] = '3-point'

    print(f"\nStart inversion: {run_name}")
    start_time = timeit.default_timer()

    if ref_ampl == 'PURA' or model_name.startswith('fixA'):
        fmin, p_out = minimizer_wrapper(0, fcn_with_progress, p_input, stat_p, meanavg=False, **kwargs)
    else:
        fmin, p_out = minimizer_wrapper(0, fcn_with_progress, p_input, stat_p, **kwargs)

    fcn_with_progress.close()
    elapsed = timeit.default_timer() - start_time
    print(f"\nElapsed time: {elapsed:.2f} seconds")
    print(f"Cost function of inversion result: {fmin.fun}")
    print('------------------')

    logname = run_name + '.txt'
    utils.writelog(runpath, logname, fmin, pars_dict, p_input, p_model, p_out)

    if save:
        with open(os.path.join(runpath, run_name + '.pkl'), 'wb') as fout:
            pickle.dump(p_out, fout, protocol=pickle.HIGHEST_PROTOCOL)

    return fmin, p_out
