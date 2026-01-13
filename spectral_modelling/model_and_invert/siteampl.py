import shutil
import os
import numpy as np
import pickle
from matplotlib import pyplot as plt
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.utils as utils
from spectral_modelling.model_and_invert.fas_log import ln_fas


def avg_mean_residual(stat_p, Z_calc):
    """
    Calculate site amplification and standard deviation in a vectorized manner,
    maintaining the original logic that filters frequencies based on a
    minimum number of events (`cfg.MIN_EVS`).
    """
    N_ev, N_sta, N_freq = stat_p.N_ev, stat_p.N_sta, stat_p.N_freq

   
    mask_3d = stat_p.M.reshape(N_ev, N_sta, N_freq)
    data_3d = stat_p.data.reshape(N_ev, N_sta, N_freq)
    model_3d = Z_calc.reshape(N_ev, N_sta, N_freq)

  
    valid_counts = np.sum(mask_3d, axis=0)  
   
    filter_mask = valid_counts >= cfg.MIN_EVS

   
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(data_3d, model_3d, out=np.ones_like(data_3d), where=mask_3d == 1)

  
    log_ratio_sum = np.sum(np.log(ratio), axis=0)
    
    
    log_gavg = np.zeros_like(valid_counts, dtype=float)
    np.divide(log_ratio_sum, valid_counts, out=log_gavg, where=filter_mask)
    
    GAVG = np.exp(log_gavg)

   
    log_ratio_sq_diff = np.square(np.log(ratio) - log_gavg)
   
    sum_sq_diff = np.sum(log_ratio_sq_diff * mask_3d, axis=0)
    
    
    gsd_variance = np.zeros_like(valid_counts, dtype=float)
    np.divide(sum_sq_diff, valid_counts, out=gsd_variance, where=filter_mask)

    GSD = np.exp(np.sqrt(gsd_variance))

    
    return GAVG.flatten(), GSD.flatten()



def fdsiteamp(runconfig, runname=cfg.RUNNAME, sptdbpath=cfg.SPTDBPATH, hcomponent='R', plot=True):
    """
    Questa funzione calcola e salva le funzioni di amplificazione di sito a(f).
    Utilizza la funzione ottimizzata avg_mean_residual.
    (Il resto della funzione rimane invariato, ma beneficia delle performance migliorate)
    """

   
    # OPEN INPUT FILEs
    model_name = utils.read_run_label(runconfig)[1]
    with open(cfg.SYNTHPARS_PATH + '/' + runname + '_' + model_name + '.pkl', 'rb') as f:
        p_model = pickle.load(f)

    N_gamma = p_model.gamma.size
    N_delta = p_model.delta.size

    subtract_noise = utils.read_run_label(runconfig)[5]
    stat_p_s, ml_s, fmin_s, fmax_s, hcompout_s = \
        utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                   hcomponent=hcomponent, component='S')
    if subtract_noise:
        stat_p_n, ml_n, fmin_n, fmax_n, hcompout_n = \
            utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                    hcomponent=hcomponent, component='N')

    N_i = stat_p_s.N_ev
    N_j = stat_p_s.N_sta
    N_k = stat_p_s.N_freq
    orids_dict = dict(zip(stat_p_s.orids, list(range(N_i))))
    stas_dict = dict(zip(stat_p_s.stas, list(range(N_j))))
    reverse_stas_dict = dict(zip(list(range(N_j)), stat_p_s.stas))

    stat_p = stat_p_s

    runpath = cfg.RUNPATH + '/' + runconfig 
    with open(runpath + '/' + runconfig + '.pkl', 'rb') as f:
        Params_slsqp = pickle.load(f)

    # extract useful parameters from Params_slsqp
    site_fi_inv = Params_slsqp.site_fi

    # build the Params object equivalent to Params_slsqp, but without eps,
    # to be used to calculate the forward model
    p_strip = utils.strip_pars(Params_slsqp, stat_p)

    ee = np.array([])
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)
    Z_calc = ln_fas(p_strip, fcn_args_ee)
    Z_calc = np.exp(Z_calc)

    if subtract_noise:
        Z_true = np.full((N_i * N_j * N_k), 1.)
        Z_diff = stat_p_s.data - stat_p_n.data
        Z_true[np.where(Z_diff > 0)] = Z_diff[np.where(Z_diff > 0)]
    else:
        Z_true = stat_p_s.data.copy()

    stat_p.set_data(Z_true)
    
    
    GAVG, GSD = avg_mean_residual(stat_p, Z_calc)


    dirpath = runpath + '/meanresiduals'
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)
    os.makedirs(dirpath)

    for s in range(len(stas_dict)):
        filename = dirpath + '/' + str(reverse_stas_dict[s]) + '_siteamp.txt'
        with open(filename, 'w') as fout:
            for k in range(N_k):
                fout.write(str(stat_p.freqs[k]) + ' ' + str(
                    GAVG[k + N_k * s]) + ' ' + str(GSD[k + N_k * s]) + '\n')

    #################
    # plot mean residuals (with variation) for each station
    # and total site amplification (site_fi * meanresidual)
    if plot is True:

        pltpath = runpath + '/plots/meanresiduals'
        if os.path.exists(pltpath):
            shutil.rmtree(pltpath)
        os.makedirs(pltpath)

        for s in range(len(stas_dict)):

            fig, ax = plt.subplots(1, 2, figsize=(8, 4), dpi=150)
            fig.suptitle(
                'sta = ' + str(reverse_stas_dict[s]) + ' , soil class = '
                + str(_c.soil_dict[reverse_stas_dict[s]]) + ', '
                + hcomponent + ', ' + runconfig, fontsize=8)

            for e in range(len(orids_dict)):
                F_e = np.zeros(N_k)
                Z_e = np.zeros(N_k)
                M_e = np.zeros(N_k)
                Z_calc_e = np.zeros(N_k)
                for k in range(N_k):
                    F_e[k] = stat_p.F[k + N_k * s + N_j * N_k * e]
                    Z_e[k] = Z_true[k + N_k * s + N_j * N_k * e]
                    M_e[k] = stat_p.M[k + N_k * s + N_j * N_k * e]
                    Z_calc_e[k] = Z_calc[k + N_k * s + N_j * N_k * e]
                index = np.where(M_e == 1)
                F_1 = F_e[index]
                Z_1 = Z_e[index]
                Z_CALC_1 = Z_calc_e[index]
                if len(F_1) > 0:
                    ax[0].loglog(F_1, (Z_1 / Z_CALC_1), c='grey', lw=0.5)
                    ax[1].loglog(F_1, (Z_1 / Z_CALC_1) * np.exp(site_fi_inv[s]),
                                 c='grey', lw=0.5)

            gavg_s = GAVG[(N_k * s): (N_k * s + N_k)]
            gsd_s = GSD[(N_k * s): (N_k * s + N_k)]
            index = np.where(gavg_s != 1.)
            ax[0].loglog(stat_p.freqs[index], gavg_s[index], c='blue', zorder=2)
            ax[0].loglog(stat_p.freqs[index], (gavg_s * gsd_s)[index], lw=0.8,
                         c='blue', zorder=2)
            ax[0].loglog(stat_p.freqs[index], (gavg_s / gsd_s)[index], lw=0.8,
                         c='blue', zorder=2)
            ax[1].loglog(stat_p.freqs[index],
                         gavg_s[index] * np.exp(site_fi_inv[s]), c='blue',
                         zorder=2)

            for i in range(2):
                ax[i].axhline(y=1., zorder=1, c='k', lw=0.8)
                ax[i].axhline(y=2.5, zorder=1, c='k', lw=0.8, ls='dotted')
                ax[i].axhline(y=0.4, zorder=1, c='k', lw=0.8, ls='dotted')
                ax[i].set_xlabel('frequency [Hz]', fontsize=7)
                ax[i].tick_params(axis='both', which='both', labelsize=7)
                ax[i].set_ylim([0.05, 40])

            ax[0].set_ylabel('residuals (FAS_real - FAS_modelled)', fontsize=7)
            ax[0].set_title('freq. dependent amplification', fontsize=7)
            ax[1].set_ylabel('residuals (FAS_real - FAS_modelled) * A',
                             fontsize=7)
            ax[1].set_title('scaled freq. dependent amplification', fontsize=7)

            plt.subplots_adjust(wspace=0.35, top=0.85, right=0.97, left=0.1,
                                bottom=0.15)
            plt.savefig(pltpath + '/' + str(reverse_stas_dict[s]) + '.png')
            plt.close(fig)
