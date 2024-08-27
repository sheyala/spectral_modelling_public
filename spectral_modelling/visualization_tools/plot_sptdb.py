import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from spectral_modelling.utils import utils
from spectral_modelling.utils import config as cfg


def plotspt(orid, sta, sptpath=cfg.SPTDBPATH, runname=cfg.RUNNAME):
    #for each component/direction, read all spectra and select the one
    #with matching evid and sta codes

    pkl_n_r = sptpath + '/' + runname +  '_noise_R.pkl'
    spt_n_r_all = utils.read_pickle_list(pkl_n_r)
    spt_n_r_stas = np.array([spt.sta for spt in spt_n_r_all])
    index_n_r_sta = np.where(spt_n_r_stas == sta)[0]
    spt_n_r_evs = np.array([spt.orid for spt in spt_n_r_all])
    index_n_r_ev = np.where(spt_n_r_evs == orid)[0]
    index = list(set(index_n_r_ev).intersection(index_n_r_sta))[0]
    spt_n_r = spt_n_r_all[index]

    pkl_n_t = sptpath + '/' + runname +  '_noise_T.pkl'
    spt_n_t_all = utils.read_pickle_list(pkl_n_t)
    spt_n_t_stas = np.array([spt.sta for spt in spt_n_t_all])
    index_n_t_sta = np.where(spt_n_t_stas == sta)[0]
    spt_n_t_evs = np.array([spt.orid for spt in spt_n_t_all])
    index_n_t_ev = np.where(spt_n_t_evs == orid)[0]
    index = list(set(index_n_t_ev).intersection(index_n_t_sta))[0]
    spt_n_t = spt_n_t_all[index]

    pkl_s_r = sptpath + '/' + runname +  '_signal_R.pkl'
    spt_s_r_all = utils.read_pickle_list(pkl_s_r)
    spt_s_r_stas = np.array([spt.sta for spt in spt_s_r_all])
    index_s_r_sta = np.where(spt_s_r_stas == sta)[0]
    spt_s_r_evs = np.array([spt.orid for spt in spt_s_r_all])
    index_s_r_ev = np.where(spt_s_r_evs == orid)[0]
    index = list(set(index_s_r_ev).intersection(index_s_r_sta))[0]
    spt_s_r = spt_s_r_all[index]

    pkl_s_t = sptpath + '/' + runname +  '_signal_T.pkl'
    spt_s_t_all = utils.read_pickle_list(pkl_s_t)
    spt_s_t_stas = np.array([spt.sta for spt in spt_s_t_all])
    index_s_t_sta = np.where(spt_s_t_stas == sta)[0]
    spt_s_t_evs = np.array([spt.orid for spt in spt_s_t_all])
    index_s_t_ev = np.where(spt_s_t_evs == orid)[0]
    index = list(set(index_s_t_ev).intersection(index_s_t_sta))[0]
    spt_s_t = spt_s_t_all[index]

    fig, axs = plt.subplots(1, 2, figsize=(10,5))
    axs[0].loglog(spt_n_r.freq, spt_n_r.amp, c='grey', alpha=0.7)
    axs[0].loglog(spt_s_r.freq, spt_s_r.amp, c='k', alpha=0.7)
    axs[0].loglog(spt_n_r.sm_freq, spt_n_r.sm_amp, c='orange', alpha=0.8)
    axs[0].loglog(spt_s_r.sm_freq, spt_s_r.sm_amp, c='r', alpha=0.8)
    axs[0].set_title('R')

    axs[1].loglog(spt_n_t.freq, spt_n_t.amp, c='grey', alpha=0.7, label='noise')
    axs[1].loglog(spt_s_t.freq, spt_s_t.amp, c='k', alpha=0.7, label='signal')
    axs[1].loglog(spt_n_t.sm_freq, spt_n_t.sm_amp, c='orange', alpha=0.8,
                  label='noise (smoothed)')
    axs[1].loglog(spt_s_t.sm_freq, spt_s_t.sm_amp, c='r', alpha=0.8,
                  label='signal (smoothed)')
    axs[1].set_title('T')

    for ax in axs:
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('FAS amplitude')

    axs[1].legend(loc='best')
    fig.suptitle('orid = ' + str(orid) + ' - sta = ' + sta)
    fig.tight_layout()
    outpath = sptpath + '/plots'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    fig.savefig(outpath + '/plot_orid' + str(orid) + '_sta' + sta + '.png')
    plt.close(fig)




# run the code for a specific evid and sta
if 0:
    myrunname = 'run1'
    myevid = 1
    mysta = 'ACOM'
    plotspt(myevid, mysta, runname=myrunname)

#run the code for all events/stations
if 0:
    myrunname = cfg.RUNNAME
    evdb = pd.read_csv(cfg.EVDB_FOR_INVERSION, sep=";")
    orids = evdb.evid_inst
    stas = evdb.stassoc
    for e in range(len(orids)):
        evstas = stas[e]
        for i in ["[", "]", "'", " "]:
            evstas = evstas.replace(i, '')
        evstas = evstas.split(',')
        for sta in evstas:
            print(orids[e], ' ', sta)
            plotspt(orids[e], sta, runname=myrunname)
