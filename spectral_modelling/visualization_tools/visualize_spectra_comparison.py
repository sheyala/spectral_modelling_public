from matplotlib import pyplot as plt
import numpy as np
import os, shutil
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils
import spectral_modelling.model_and_invert.fas_log as fas
import pickle

soil_dict = _c.soil_dict

#######################################
# This program creates plots of the spectra (original, smoothed, inverted
# without epsilons, inverted with epsilons - with associated uncertainty)

def plot_spectra_comparison(runname=cfg.RUNNAME, 
                            sptdbpath=cfg.SPTDBPATH, runpath=cfg.RUNPATH,
                            hcomponent='R',
                            model_name='malagniniQ0', method='SLSQP', 
                            use_uncert='noeps', ref_ampl='meanA', 
                            subtract_noise=False, weights=None, 
                            bounds='NE_Italy'):
    """
    Create comparison plots with the recorded spectrum and the modelled one
    for each available station and event
    """

    run_name = utils.set_run_label(model_name, use_uncert, weights, bounds,
                                    subtract_noise, ref_ampl, method,
                                    hcomponent)


    ########################################
    # OPEN INPUT FILE
    stat_p, ml, mag, fmin, fmax = \
        utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                   hcomponent=hcomponent)

    stat_p_n, ml_n, mag_n, fmin_n, fmax_n =  \
        utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                   hcomponent=hcomponent, component='N')

    Z_signal = stat_p.data
    Z_noise = stat_p_n.data

    orids = stat_p.orids
    stas = stat_p.stas
    freqs = stat_p.freqs
    N_i = stat_p.N_ev
    N_j = stat_p.N_sta
    N_k = stat_p.N_freq
    orids_dict = dict(zip(orids,list(range(N_i))))
    stas_dict = dict(zip(stas,list(range(N_j))))
    ml_dict = dict(zip(orids,ml))
    reverse_orids_dict = dict(zip(list(range(N_i)),orids))
    reverse_stas_dict = dict(zip(list(range(N_j)),stas))


    dirpath = cfg.SYNTHPARS_PATH
    with open(dirpath + '/' + runname + '_' + model_name + '.pkl', 'rb') as f:
        p_true = pickle.load(f)
    N_gamma = p_true.gamma.size
    N_delta = p_true.delta.size

    with open(runpath + '/' + run_name + '/' + run_name +'.txt', 'r') as fname:
        pars_slsqp = []
        for line in fname:
            if not line.startswith("#"):
                line = line.strip()
                line = line.split()
                pars_slsqp.append(float(line[3]))

    pars_slsqp = np.array(pars_slsqp)

    #extract useful parameters from pars_slsqp
    alpha_inv = pars_slsqp[0:N_i]
    mag_inv = utils.m0_to_mag_HK(np.exp(alpha_inv)) #Mw
    ml_inv = utils.mag_to_ml_munafo(mag_inv) #ML
    

    #build the Params object equivalent to pars_slsqp, but without eps, 
    #to be used to calculate the forward model
    pars_slsqp_Pars = utils.array_to_Param(pars_slsqp, stat_p, N_gamma, N_delta)
    p_strip = utils.strip_pars(pars_slsqp_Pars, stat_p)

    ee = np.array([])
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)
    Z_calc = fas.ln_fas(p_strip, fcn_args_ee)
    Z_calc = np.exp(Z_calc)

    Z_true = Z_signal.copy()
    if subtract_noise == True:
        Z_true = np.full((N_i*N_j*N_k),1.)
        Z_diff = Z_signal - Z_noise
        Z_true[np.where(Z_diff > 0)] = Z_diff[np.where(Z_diff > 0)]

    plotpath = runpath + '/' + run_name  + '/plots'
    if (os.path.exists(plotpath) == 0):
        os.makedirs(plotpath)

    ##############################
    # plot models output of SLSQP inversion
    if 1:
        outpath = plotpath + '/spectra_' + run_name
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.makedirs(outpath)

        for e in range(len(orids_dict)):
            for s in range(len(stas_dict)):
                F_DST = np.zeros(N_k)
                Z_DST = np.zeros(N_k)
                Z_calc_DST = np.zeros(N_k)
                M_DST = np.zeros(N_k)
                for k in range(N_k):
                    F_DST[k] = stat_p.F[k+N_k*s+N_j*N_k*e]
                    Z_DST[k] = Z_true[k+N_k*s+N_j*N_k*e]
                    Z_calc_DST[k] = Z_calc[k+N_k*s+N_j*N_k*e]
                    M_DST[k] = stat_p.M[k+N_k*s+N_j*N_k*e]
                F_DST_1 = F_DST[np.where(M_DST==1)]
                Z_DST_1 = Z_DST[np.where(M_DST==1)]
                Z_calc_DST_1 = Z_calc_DST[np.where(M_DST==1)]

                if len(F_DST_1) > 0:
                    orid = reverse_orids_dict[e]
                    sta = reverse_stas_dict[s]

                    pkl_s_r = sptdbpath + '/' + runname +  '_signal_R.pkl'
                    spt_s_r_all = utils.read_pickle_list(pkl_s_r)
                    spt_s_r_stas = np.array([spt.sta for spt in spt_s_r_all])
                    index_s_r_sta = np.where(spt_s_r_stas == sta)[0]
                    spt_s_r_evs = np.array([spt.orid for spt in spt_s_r_all])
                    index_s_r_ev = np.where(spt_s_r_evs == orid)[0]
                    index = list(set(index_s_r_ev).intersection(index_s_r_sta))[0]
                    spt_s_r = spt_s_r_all[index]

                    index = np.where(spt_s_r.freq >= 0.5)
                    F_db_v0 = spt_s_r.freq[index]
                    Z_db_v0 = spt_s_r.amp[index]
                    index = np.where(spt_s_r.freq <= 25)
                    F_db_v1 = F_db_v0[index]
                    Z_db_v1 = Z_db_v0[index]

                    plt.xscale('log',base=10.)
                    plt.yscale('log',base=10.)
                    plt.plot(F_db_v1, Z_db_v1, c='dimgrey', zorder=-1, lw=1)
                    plt.plot(F_DST_1, Z_DST_1, label='Z_true', c = 'blue')
                    plt.plot(F_DST_1, Z_calc_DST_1, c='green', label = 'Inverted FAS', zorder=1, lw=2)

                    plt.title('ev = %s; sta = %s; soil class =%s; M$_L$ = %.1f; R$_{hyp}$ = %.0f km' % (reverse_orids_dict[e], reverse_stas_dict[s], soil_dict[reverse_stas_dict[s]], ml_inv[e], (stat_p.R[N_j*N_k*e + N_k*s]/100000)))
                    plt.legend()
                    plt.savefig(outpath+'/orid'+str(reverse_orids_dict[e])+'_'+str(reverse_stas_dict[s])+'.png')
                    plt.close()



##############################
