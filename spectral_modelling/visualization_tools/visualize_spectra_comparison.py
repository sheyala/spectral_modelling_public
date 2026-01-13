from matplotlib import pyplot as plt
import numpy as np
import os, shutil
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils
import spectral_modelling.model_and_invert.fas_log as fas
import pickle

soil_dict = _c.soil_dict

def plot_spectra_comparison(runname=cfg.RUNNAME,
                            sptdbpath=cfg.SPTDBPATH, runpath=cfg.RUNPATH,
                            hcomponent='R',
                            model_name='malagniniQ0', method='SLSQP',
                            use_uncert='noeps', ref_ampl='meanA',
                            subtract_noise=False, weights=None,
                            bounds='NE_Italy'):
 

    run_name = utils.set_run_label(model_name, use_uncert, weights, bounds,
                                     subtract_noise, ref_ampl, method,
                                     hcomponent)

  
    stat_p, ml, mag, fmin, fmax = \
        utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                     hcomponent=hcomponent)

 
    results_path = os.path.join(runpath, run_name, f"{run_name}.pkl")
    with open(results_path, 'rb') as f:
        p_out = pickle.load(f)

    N_i = stat_p.N_ev
    N_j = stat_p.N_sta
    N_k = stat_p.N_freq
    
    reverse_orids_dict = dict(zip(range(N_i), stat_p.orids))
    reverse_stas_dict = dict(zip(range(N_j), stat_p.stas))

   
    ee = np.array([])
    N_gamma = p_out.gamma.size
    N_delta = p_out.delta.size
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)
    Z_calc = fas.ln_fas(p_out, fcn_args_ee)
    Z_calc = np.exp(Z_calc)

   
    Z_true = stat_p.data.copy()
    if subtract_noise:
        stat_p_n, _, _, _, _ = \
            utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                       hcomponent=hcomponent, component='N')
        Z_diff = stat_p.data - stat_p_n.data
       
        Z_true = np.where(Z_diff > 0, Z_diff, 1.0)

   
    print("Loading spectra database")
    pkl_s_r = os.path.join(sptdbpath, f"{runname}_signal_R.pkl")
    spt_s_r_all = utils.read_pickle_list(pkl_s_r)
    
    spt_dict = {(spt.orid, spt.sta): spt for spt in spt_s_r_all}
    # -------------------------------------------------------------------------

    outpath = os.path.join(runpath, run_name, 'plots', f'spectra_{run_name}')

    
    if os.name == 'nt':
        outpath = os.path.abspath(outpath)
        if not outpath.startswith("\\\\?\\"):
            outpath = "\\\\?\\" + outpath
    # ----------------------------------

    os.makedirs(outpath, exist_ok=True)

    print("Start graph creation...")
    for e in range(N_i):
        for s in range(N_j):
            orid = reverse_orids_dict[e]
            sta = reverse_stas_dict[s]

            output_filename = os.path.join(outpath, f'orid{orid}_{sta}.png')
            
            if os.path.exists(output_filename):
                continue

            
            idx_start = e * N_j * N_k + s * N_k
            idx_end = idx_start + N_k

            M_DST = stat_p.M[idx_start:idx_end]
            
            if np.sum(M_DST) == 0:
                continue

            
            F_DST_1 = stat_p.F[idx_start:idx_end][M_DST==1]
            Z_DST_1 = Z_true[idx_start:idx_end][M_DST==1]
            Z_calc_DST_1 = Z_calc[idx_start:idx_end][M_DST==1]

           
            spt_s_r = spt_dict.get((orid, sta))
            if spt_s_r is None:
               
                continue
            # ----------------------------------------------------------------------

            
            index = (spt_s_r.freq >= 0.5) & (spt_s_r.freq <= 25)
            F_db_v1 = spt_s_r.freq[index]
            Z_db_v1 = spt_s_r.amp[index]

            # --- PLOTTING
            plt.figure(figsize=(8, 6)) 
            plt.xscale('log', base=10.)
            plt.yscale('log', base=10.)

            
            plt.plot(F_db_v1, Z_db_v1, c='dimgrey', zorder=-1, lw=1)

            
            plt.plot(F_DST_1, Z_DST_1, label='Z_true', c = 'blue')

            
            plt.plot(F_DST_1, Z_calc_DST_1, c='green', label = 'Inverted FAS', zorder=1, lw=2)

            
            current_soil = soil_dict.get(sta, "N/D") 
            current_ml = ml[e]
            
            current_rhyp = stat_p.R[idx_start] / 100000

           
            plt.title('ev = %s; sta = %s; soil class =%s; M$_L$ = %.1f; R$_{hyp}$ = %.0f km' % (
                orid, sta, current_soil, current_ml, current_rhyp))

            plt.legend()
            plt.savefig(output_filename, bbox_inches='tight')
            plt.close()

    print("Graph creation completed.")


##############################
