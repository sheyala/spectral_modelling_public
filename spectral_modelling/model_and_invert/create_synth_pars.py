import numpy as np
import pandas as pd
import pickle
import os
import spectral_modelling.utils.utils as utils
import spectral_modelling.utils.myClass as myC
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.constants as _c

def synthetic_p_in(runname=cfg.RUNNAME, sptdbpath=cfg.SPTDBPATH,
                   hcomponent='R', model='malagniniQ0',
                   evdb=cfg.EVDB_FOR_INVERSION, stadb=cfg.STADB_FOR_INVERSION, 
                   assocdb=cfg.ASSOCDB):
    """
    Create sample input parameter set(s) and store it as pickle object
    """

    # OPEN INPUT FILE and read information

    stat_p, ml, mag, fmin, fmax = \
        utils.read_input_component(runname=runname, sptdbpath=sptdbpath,
                                   hcomponent=hcomponent)
    N_i = stat_p.N_ev #n. of events in database
    N_j = stat_p.N_sta #n. of stations in database

    m0_calc = utils.mag_to_m0_HK(mag)
    alpha_calc = np.log(m0_calc)
    beta_calc = utils.m0_to_beta_Brune(m0_calc, _c.AVG_STDROP)


    # parse available associations, events and stations
    assocdb = pd.read_csv(assocdb, sep=";")
    evdb = pd.read_csv(evdb, sep=";")
    stadb = pd.read_csv(stadb, sep=";")

    evids_inst = evdb.evid_inst.unique()
    stas = stadb.sta.unique()
    evids_dict = dict(zip(evids_inst, list(range(N_i))))
    stas_dict = dict(zip(stas, list(range(N_j))))

    dirpath = cfg.SYNTHPARS_PATH
    if os.path.exists(dirpath) == 0:
        os.makedirs(dirpath)

    # build and save a set of model parameters
    # "malagnini_Q0": fc from stress drop; deltas from Malagnini;
    #                 delta[1] fixed to 0
    match model:

        case 'malagniniQ0':
            alpha = alpha_calc
            beta = beta_calc
            gamma = np.array([-1., -1.6, -1.2, -1.3, -0.5, -0.95, -1.2, -1.8, 
                              -1.2, -0.5], dtype='d')
            delta = np.array([500., 0.], dtype='d')
            site_fi = np.full(N_j, 0.)
            site_k = np.full(N_j, 0.04)
            eps_source = np.zeros(N_i)
            eps_path = np.array([0.])
            eps_site = np.zeros(N_j)
            synthpars = myC.Params()
            synthpars.set_pars(alpha, beta, gamma, delta, site_fi, site_k, 
                               eps_source, eps_path, eps_site)

            N_gamma = gamma.size
            synthpars.vary_all(True)
            for i in range(N_gamma):
                synthpars.set_vary((2 * N_i + i,), False)
            synthpars.set_vary((2 * N_i + N_gamma + 1,), False)

        case _:
            raise ValueError(f"Model {model} not implemented")

    with open(dirpath + '/' + runname + '_' + model + '.pkl', 'wb') as out:
        pickle.dump(synthpars, out, protocol=pickle.HIGHEST_PROTOCOL)
