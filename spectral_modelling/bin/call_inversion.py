import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils
from spectral_modelling.model_and_invert.create_synth_pars import synthetic_p_in
from spectral_modelling.model_and_invert.inversion import fas_invert
from spectral_modelling.model_and_invert.siteampl import fdsiteamp
from spectral_modelling.visualization_tools.visualize_spectra_comparison import plot_spectra_comparison


def call_inv_modules(runname=cfg.RUNNAME, sptdbpath=cfg.SPTDBPATH, 
                     outpath=cfg.RUNPATH, hcomponent='R', 
                     model_name='malagniniQ0', method='SLSQP', jac='exact', 
                     use_uncert='noeps', ref_ampl='meanA', 
                     subtract_noise=False, weights=None, bounds='NE_Italy', 
                     save=True, calc_synth=False, plots=False):

    if calc_synth:
        synthetic_p_in()

    fmin_noe, p_out_noe = fas_invert(runname=runname, sptdbpath=sptdbpath, 
                                     outpath=outpath, hcomponent=hcomponent, 
                                     model_name=model_name, method=method, 
                                     jac=jac, use_uncert=use_uncert, 
                                     ref_ampl=ref_ampl,  
                                     subtract_noise=subtract_noise,
                                     weights=weights, bounds=bounds, save=save)

    configname = utils.set_run_label(model_name, use_uncert, weights, bounds, 
                                   subtract_noise, ref_ampl, method, hcomponent)
    print(configname)

    if plots:
        plot_spectra_comparison(runname=runname, sptdbpath=sptdbpath, 
                                runpath=outpath,hcomponent=hcomponent, 
                                model_name=model_name, method=method,
                                use_uncert=use_uncert, ref_ampl=ref_ampl,  
                                subtract_noise=subtract_noise,
                                weights=weights, bounds=bounds)


    fdsiteamp(configname, runname=runname, sptdbpath=sptdbpath,
              hcomponent=hcomponent, plot=True)



#run the inversion using default configuration, except that it calculates 
#the initial synthetic parameters and plots the result
if 1:
    call_inv_modules(calc_synth=True, plots=True)

