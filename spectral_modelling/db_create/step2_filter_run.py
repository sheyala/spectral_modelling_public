import pandas as pd
import numpy as np
from collections import OrderedDict
from collections import Counter
import os
import pickle
from spectral_modelling.utils import config as cfg

def filter_evstas(inpath_wf=cfg.WFDBPATH,outpath=cfg.DBPATH,dbname=cfg.DBNAME,
                  eventdb=cfg.EVDB, stationdb=cfg.STADB, runname=cfg.RUNNAME, 
                  snrmin=cfg.SNRMIN, nevsmin=cfg.MIN_EVS, nstamin=cfg.MIN_STAS, 
                  min_hypodist=cfg.MIN_HYPODIST, max_hypodist=cfg.MAX_HYPODIST,
                  verbose=0):
    """
    This script filters the input lists of events and stations by keeping only 
    stations with at least nevsmin records and events recorded by at least 
    nstamin stations. 
    It also apply a filter on the signal-to-noise ratio level and on the 
    hypocentral distance range.

    The input lists of available events and stations must be in the form
    (EVID;ORID;ETIME;ELAT;ELON;EDEPTH;ML) and (STA;SLAT;SLON;SELEV) 
    respectively, e.g.:

    evid;orid;etime;elat;elon;edepth;ml
    9;20134158;2013-10-07T17:35:52.544Z;63.8574;-22.4215;4.92578;2.9  

    sta;slat;slon;selev
    ISS;63.86293;-22.30294;71.0   
    
    inpath_wf  path to input waveforms
    outpath    path to store results
    eventdb    csv file containing info on inital pool of events
    stationdb  csv file containing info on inital pool of stations
    runname    name of the run (associated to choices on constraints)
    snrmin     min signal-to-noise ratio for each record
    nevsmin    min n. of events recorded by each station
    nstamin    min n. of stations for each event
    min_hypodist, max_hypodist minimum and maximum hypocentral distance to be 
    used
    """

    #------------------------------------------------------------------
    # PREPARE OUTPUT 
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    f_out_ev = os.path.join(outpath, dbname + "_events_" + runname + ".txt")
    f_out_sta = os.path.join(outpath, dbname + "_stas_" + runname + ".txt")
    f_out_readme =os.path.join(outpath, "README_" + dbname + "_" + runname 
                               + ".txt")

    #------------------------------------------------------------------
    # PARSE INPUT DATA AND PARAMETERS

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
    # extract evid and corresponding evid_inst list
    evids = evdb.evid
    evids_inst = evdb.evid_inst

    #------------------------------------------------------------------
    # BUILD INITIAL POOL OF MATCHES (evid_inst+sta), STAS, EVS
    init_match = []
    init_sta = []
    init_ev = []
    for e in range(len(evids)):
        try:
            with open(inpath_wf + '/' + str(evids_inst[e]) + '.pkl', 'rb') as f:
                wfpkl = pickle.load(f)
            if verbose:
                print(len(wfpkl['stacode']), ' stas found for evid ', 
                      str(evids_inst[e]))
        except:
            if verbose:
                print("event with evid ", evids_inst[e], "not in inpath_wf")
            continue

        # index to filter by snr score
        indR = np.where(np.array(wfpkl['snr']['R'], copy=False) > snrmin)[0]
        indT = np.where(np.array(wfpkl['snr']['T'], copy=False) > snrmin)[0]
        ind_snr = list(set(indT) & set(indR))

        # index to filter by hypocentral distance
        indD1 = np.where(np.array(wfpkl['hypodist'], copy=False) > min_hypodist)[0]
        indD2 = np.where(np.array(wfpkl['hypodist'], copy=False) < max_hypodist)[0]
        ind_hypo = list(set(indD1) & set(indD2))

        ind_filt = list(set(ind_snr) & set(ind_hypo))
        filt_sta = [wfpkl['stacode'][i] for i in ind_filt]

        for j in range(len(filt_sta)):
            init_match.append([evids[e],filt_sta[j]])
            init_sta.append(filt_sta[j])
            init_ev.append(evids[e])

    #------------------------------------------------------------------
    # RECURSIVE FILTER
    # recursive filter on n. of events recorded by each station (nevsmin)
    # and n. of recordings for each event (nstamin)
    prev_match = init_match
    prev_sta = init_sta
    prev_ev = init_ev
    prev_filtmatch = {}
    j = 0
    updated_status = False
    while(j<10):
        #------------------------------------------------------------------
        # filter on origins
        ev_count = Counter(list(prev_ev))
        filtmatch = {}
        for i in ev_count: 
            if(ev_count[i]>=nstamin):
                filtmatch[i] = ev_count[i]
        # check if filtering is finished and, if so, exit cycle
        if(filtmatch == prev_filtmatch):
            break
        temp_match_1 = []
        temp_sta_1 = []
        temp_ev_1 = []
        for i in range(len(prev_match)):
            if(prev_ev[i] in filtmatch):
                temp_match_1.append(prev_match[i])
                temp_sta_1.append(prev_match[i][1])
                temp_ev_1.append(prev_match[i][0])
        #------------------------------------------------------------------
        # filter on stations
        sta_count = Counter(list(temp_sta_1))
        filtsta = {}
        for i in sta_count: 
            if(sta_count[i]>=nevsmin):
                 filtsta[i] = sta_count[i]
        temp_match_2 = []
        temp_sta_2 = []
        temp_ev_2 = []
        for i in range(len(temp_match_1)):
            if(temp_sta_1[i] in filtsta):
                temp_match_2.append(temp_match_1[i])
                temp_sta_2.append(temp_match_1[i][1])
                temp_ev_2.append(temp_match_1[i][0])
        #------------------------------------------------------------------
        # store outputs for next iteration
        prev_match = temp_match_2
        prev_sta = temp_sta_2
        prev_ev = temp_ev_2
        prev_filtmatch = filtmatch
        j = j+1

    #------------------------------------------------------------------
    # SAVE OUTPUT
    print('n. of iterations ', j)
    if (filtmatch != {}):
        updated_status = True

        # create dicts of events and stations from filt match
        matchevs =  {id:[] for id in prev_ev}
        for pmi in prev_match:    
            matchevs[pmi[0]].append(pmi[1])
        matchstas =  {id:[] for id in prev_sta}
        for pmi in prev_match:    
            matchstas[pmi[1]].append(pmi[0])

        def_ev = list(OrderedDict.fromkeys(prev_ev))
        df_ev = pd.DataFrame(columns=["evid","evid_inst","etime","elat",
                                      "elon", "edepth", "ml", "mag", 
                                      "nstassoc", "stassoc"])
        for i in range(len(def_ev)):
            current_ev = evinv[def_ev[i]]
            df_ev.loc[i] = pd.Series({"evid":def_ev[i],
                                      "evid_inst":str(current_ev[0]),
                                      "etime":str(current_ev[1]),
                                      "elat":str(current_ev[2]),
                                      "elon":str(current_ev[3]),
                                      "edepth":str(current_ev[4]),
                                      "ml":str(current_ev[5]),
                                      "mag":str(current_ev[6]),
                                      "nstassoc":len(matchevs[def_ev[i]]),
                                      "stassoc":matchevs[def_ev[i]]})
        df_ev.to_csv(f_out_ev, sep=";", index=False)


        def_sta = list(OrderedDict.fromkeys(prev_sta))
        df_sta = pd.DataFrame(columns=["net","sta","slat","slon","nevassoc"])
        for i in range(len(def_sta)):
            current_sta = stainv[def_sta[i]]
            df_sta.loc[i] = pd.Series({"net":str(current_sta[0]),
                            "sta":str(current_sta[1]),
                            "slat":str(current_sta[2]),
                            "slon":str(current_sta[3]),
                            "nevassoc":len(matchstas[def_sta[i]])})
        df_sta.to_csv(f_out_sta, sep=";", index=False)

        with open(f_out_readme,"w") as f:
            f.write(f'{"#README: run name = ":<20}{runname}\n')
            f.write(f'{"#snrmin":<9}{"min_hypodist":<14}{"max_hypodist":<14}{"nstamin":<9}{"nevsmin":<9}\n')
            f.write(f'{snrmin:<9}{min_hypodist:<14}{max_hypodist:<14}{nstamin:<9}{nevsmin:<9}\n')

    return updated_status

