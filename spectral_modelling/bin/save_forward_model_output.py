import numpy as np
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.myClass as myC
import spectral_modelling.utils.utils as utils
import spectral_modelling.model_and_invert.fas_log as fas

gamma_malagnini = np.array([-1., -1.6, -1.2, -1.3, -0.5, -0.95, -1.2, -1.8, 
                            -1.2, -0.5], dtype='d')
def single_spectrum_pars(alpha, beta, delta, site_fi, site_k, 
                         eps_source, eps_path, eps_site, gamma=gamma_malagnini):
#   alpha = log(M0 [dyn * cm])
#   beta = fc [Hz]
#   delta = ([Q0, eta_Q])
#   site_fi = freq.-independent site scaling
#   site_k = k-value at site [s]
#   eps_source = uncertainty on ln(source part)
#   eps_path = uncertainty on ln(path part)
#   eps_site = uncertainty on ln(site part
#   gamma = geom. spreading function coefficients
#           [by default: from Malagnini et al.])
    
    scalarcondition = np.isscalar(alpha), np.isscalar(beta), \
                      np.isscalar(site_fi), np.isscalar(site_k), \
                      np.isscalar(eps_source), np.isscalar(eps_path), \
                      np.isscalar(eps_site)
    if False in scalarcondition:
        print('ERROR All parameters (except delta and gamma) must be scalar values')
        return None

    pars = myC.Params()
    pars.set_pars(np.array([alpha]), np.array([beta]), gamma, delta, 
                  np.array([site_fi]), np.array([site_k]), 
                  np.array([eps_source]), np.array([eps_path]), 
                  np.array([eps_site]))
    return pars


def save_fwmodel(spt_pars, hypodist, outfile, staname=None):
    """
    This program saves the output of the forward model 
    FOR AN INDIVIDUAL SPECTRUM in csv format
    
    spt_pars = myClass.Params() instance for one event, one station
    hypodist = hypocentral distance [km]
    outfile = absolute path to save output
    staname = station name (if you want to add the site amplification term)
    """

    N_gamma = spt_pars.gamma.size
    N_delta = spt_pars.delta.size
    
    #build the Params object equivalent to spt_pars, but without eps, 
    #to be used to calculate the forward model
    F = cfg.FREQS
    R = np.full_like(F, hypodist*1.e+5)
    stat_p = myC.StaticParams(F=F, R=R, N_ev=1, N_sta=1, N_freq=len(F))
    p_strip = utils.strip_pars(spt_pars, stat_p)

    ee = np.array([])
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)
    Z_calc = fas.ln_fas(p_strip, fcn_args_ee)
    Z_calc_rock = np.exp(Z_calc)
    Z_calc_site = np.exp(Z_calc)

    #if station has associated site response, add it to the forward model
    #for the calculation of Z_calc_site
    if staname:
        try:
            stat_p_all, ml, mag, fmin, fmax = \
                utils.read_input_component(runname=cfg.RUNNAME, 
                                        sptdbpath=cfg.SPTDBPATH,
                                        hcomponent='R')
            stas = stat_p_all.stas
            N_j = stat_p_all.N_sta
            stas_dict = dict(zip(stas,list(range(N_j))))

            if staname in stas_dict:
                run_name = utils.set_run_label(model_name='malagniniQ0', 
                                            use_uncert='noeps', weights=None, 
                                            bounds='NE_Italy', subtract_noise=False, 
                                            ref_ampl='meanA', method='SLSQP',
                                            hcomponent='R')
                fdamplpath = cfg.RUNPATH + '/' + run_name + '/meanresiduals'
                fname = staname + '_siteamp.txt'
                fdampl = np.loadtxt(fdamplpath+'/'+fname, usecols=(1,))
                Z_calc_site = Z_calc + np.log(fdampl)
                Z_calc_site = np.exp(Z_calc_site)
            else:
                print('Site amplification for station ' + staname + 
                    ' not found - FAS model for rock site will be used \n')
        except:
            print('Site amplifications related to run "' + 
                  run_name + 
                  '" not found - FAS model for rock site will be used \n')


    with open(outfile, 'w+') as fout:
        fout.write('#freq[Hz],FAS_calc_rock[cm],FAS_calc_site[cm]\n')
        for i in range(len(F)):
            fout.write(str(F[i]) + ',' + 
                       str(Z_calc_rock[i]) + ',' +
                       str(Z_calc_site[i]) + '\n')




##############################
# EXAMPLE
# NOTE:
# spt_pars = single_spectrum_pars(alpha = log(M0 [dyn * cm]), 
#                                 beta = fc [Hz], 
#                                 delta = ([Q0, eta_Q]), 
#                                 site_fi = freq.-independent site scaling, 
#                                 site_k = k-value at site [s], 
#                                 eps_source = uncertainty on ln(source part), 
#                                 eps_path = uncertainty on ln(path part), 
#                                 eps_site = uncertainty on ln(site part),
#                                 gamma = geom. spreading function coefficients
#                                         [by default: from Malagnini et al.])
#
# output velocity FAS is in cm
#
# NOTE the run_name from where to read site amplification is hardcoded in
# save_fwmodel() and has to be modified inside the function, if needed
            
if 0:
    outpath = '/home/laura/Desktop/per_shakelab/spectral_modelling_public/spectral_modelling/data/OGS2019/run1_output'
    outfile = outpath + '/ML3.82_HD19_GEPF.txt'

    spt_pars = single_spectrum_pars(49.841, 1.814, np.array([500.,0.]), 
                                    0.059, 0.038, 0., 0., 0.)
    hypodist=19. # in km
    staname = 'GEPF'
    save_fwmodel(spt_pars, hypodist, outfile, staname=staname)




