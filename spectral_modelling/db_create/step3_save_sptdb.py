import pandas as pd
from obspy import Stream
import os
import pickle
from spectral_modelling.utils import config as cfg
from spectral_modelling.utils import utils
from spectral_modelling.utils.myClass import Spectrum


def create_sptpkl(inpath_wf=cfg.WFDBPATH, outpath=cfg.SPTDBPATH, 
                  runname=cfg.RUNNAME, eventdb=cfg.EVDB_FOR_INVERSION,
                  stationdb=cfg.STADB_FOR_INVERSION, freqlims=cfg.FREQLIMS):
    """
    Read filtered lists of events and associated stations, extract the corresponding processed waveforms and store their FAS (Fourier amplitude spectra) in the form of pickle files containing Spectrum() objects for radial/transverse components and for both signal and noise. 

    inpath_wf  path to input waveforms
    outpath    path to store results
    eventdb    txt file containing info on filtered pool of events
    stationdb  txt file containing info on filtered pool of stations
    runname    name of the run (associated to choices on constraints)
    """

    #------------------------------------------------------------------
    # PARSE INPUT EVENTS AND STATIONS
    # parse available events and stations
    evdb = pd.read_csv(eventdb, sep=";")
    stadb = pd.read_csv(stationdb, sep=";")
    # create corresponding dictionaries
    evinv = {}
    for j in range(len(evdb)):
        evinv[evdb.iloc[j].evid] = evdb.iloc[j].values[1:].tolist()
    stainv = {}
    for j in range(len(stadb)):
        stainv[stadb.iloc[j].sta] = stadb.iloc[j].values.tolist()
    # extract evid and corresponding evids_inst list
    evids = evdb.evid
    evids_inst = evdb.evid_inst

    #------------------------------------------------------------------
    # EXTRACT CORRESPONDING WF
    # store all data in four Streams (signal/noise and R/T component)
    wf_signal_r = Stream()
    wf_signal_t = Stream()
    wf_noise_r = Stream()
    wf_noise_t = Stream()    
    wf_evcodes = []
    wf_orcodes = []
    wf_hypo = []

    for e in range(len(evids)):
        try:
            with open(inpath_wf + '/' + str(evids_inst[e]) + '.pkl', 'rb') as f:
                wfpkl = pickle.load(f)        
        except:
            print("event with evid ", evids_inst[e], "not in inpath_wf")
            continue

        signal_r = wfpkl['signal']['R']
        signal_t = wfpkl['signal']['T']
        noise_r = wfpkl['noise']['R']
        noise_t = wfpkl['noise']['T']
        stacode = wfpkl['stacode']
        hypodist = wfpkl['hypodist']

        assoc_stas = evinv[evids[e]][-1]
        for i in ["[", "]", "'", " "]:
            assoc_stas = assoc_stas.replace(i, '')
        assoc_stas = assoc_stas.split(',')

        # indexes of stations appearing in the associated list for current event
        index = [i for i, s in enumerate(stacode) if s in assoc_stas]
        for i in index:
            wf_signal_r.append(signal_r[i])
            wf_signal_t.append(signal_t[i])
            wf_noise_r.append(noise_r[i])
            wf_noise_t.append(noise_t[i])
            wf_evcodes.append(evids[e])
            wf_orcodes.append(evids_inst[e])
            wf_hypo.append(hypodist[i])


    #------------------------------------------------------------------
    # CALCULATE SPECTRA AND STORE IN PICKLE OUTPUT
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    out_signal_r = open(outpath + '/' + runname + '_signal_R.pkl', 'wb')
    out_signal_t = open(outpath + '/' + runname + '_signal_T.pkl', 'wb')
    out_noise_r = open(outpath + '/' + runname + '_noise_R.pkl', 'wb')
    out_noise_t = open(outpath + '/' + runname + '_noise_T.pkl', 'wb')

    for i in range(len(wf_evcodes)):

        evid = wf_evcodes[i]
        orid = wf_orcodes[i]
        sta = wf_signal_r[i].stats['station']
        chan = wf_signal_r[i].stats['channel'][-1]
        hypdist = wf_hypo[i]
        ml = evinv[evid][5]
        mag = evinv[evid][6]

        # ----------------------------------------------------------------------
        # for signal, R
        spt_signal_r = Spectrum(orid=orid, sta=sta, chan=chan, 
                                hypdist=hypdist, ml=ml, mag=mag,
                                freq_lims=freqlims)
        utils.fas_calc_and_store(wf_signal_r, spt_signal_r, i)
        pickle.dump(spt_signal_r, out_signal_r, 
                    protocol=pickle.HIGHEST_PROTOCOL)

        # ----------------------------------------------------------------------
        # for signal, T
        spt_signal_t = Spectrum(orid=orid, sta=sta, chan=chan, 
                                hypdist=hypdist, ml=ml, mag=mag, 
                                freq_lims=freqlims)
        utils.fas_calc_and_store(wf_signal_t, spt_signal_t, i)
        pickle.dump(spt_signal_t, out_signal_t, 
                    protocol=pickle.HIGHEST_PROTOCOL)

        # ----------------------------------------------------------------------
        # for noise, R
        spt_noise_r = Spectrum(orid=orid, sta=sta, chan=chan,
                               hypdist=hypdist, ml=ml, mag=mag, 
                               freq_lims=freqlims)
        utils.fas_calc_and_store(wf_noise_r, spt_noise_r, i)
        pickle.dump(spt_noise_r, out_noise_r, 
                    protocol=pickle.HIGHEST_PROTOCOL)

        # ----------------------------------------------------------------------
        # for noise, T
        spt_noise_t = Spectrum(orid=orid, sta=sta, chan=chan, 
                               hypdist=hypdist, ml=ml, mag=mag, 
                               freq_lims=freqlims)
        utils.fas_calc_and_store(wf_noise_t, spt_noise_t, i)
        pickle.dump(spt_noise_t, out_noise_t, 
                    protocol=pickle.HIGHEST_PROTOCOL)

    out_signal_r.close()
    out_signal_t.close()
    out_noise_r.close()
    out_noise_t.close()

