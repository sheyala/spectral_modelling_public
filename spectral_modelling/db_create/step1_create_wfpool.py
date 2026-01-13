import pickle
import pandas as pd
import os
from obspy import Stream, read, read_inventory
from obspy.core.utcdatetime import UTCDateTime
from spectral_modelling.utils import config as cfg
from spectral_modelling.utils import utils
import matplotlib.pyplot as plt

#TODO
#make sure that outcomes are always in velocity (in remove_response)


def create_wfpool(dataformat, inpath_wf=cfg.WFPATH, outpath=cfg.WFDBPATH, 
                  evdb=cfg.EVDB, stadb=cfg.STADB, assocdb=cfg.ASSOCDB, 
                  tdura=cfg.TDURA, removeresp=True, respfiles=cfg.RESPPATH,
                  verbose=0):
    match dataformat:
        case 'sac' | 'SAC':
            create_wfpool_sac(inpath_wf=inpath_wf, outpath=outpath, 
                              evdb=evdb, stadb=stadb, assocdb=assocdb, 
                              tdura=tdura, verbose=verbose)
        case 'mseed' | 'MSEED' | 'miniseed':
            create_wfpool_mseed(inpath_wf=inpath_wf, outpath=outpath, 
                                evdb=evdb, stadb=stadb, assocdb=assocdb, 
                                tdura=tdura, removeresp=removeresp, 
                                respfiles=respfiles, verbose=verbose)
        case _:
            print('Unknown data format - only sac or mseed allowed')



def create_wfpool_sac(inpath_wf=cfg.WFPATH, outpath=cfg.WFDBPATH, evdb=cfg.EVDB,
                      stadb=cfg.STADB, assocdb=cfg.ASSOCDB, tdura=cfg.TDURA, verbose=0):
    """
    Read processed waveforms for selected events and stations 
    and store information separately for each event in the form of pickle files 
    containing Stream() objects with Radial/Transverse components,
    for both signal and data.
    Duration of signal and noise windows is currently equal and chosen through
    the tdura parameter, as (tP-tdura-1, tP-1) for noise and (tS, tS+tdura) 
    for signal)
    """

    # parse available associations, events and stations
    assocdb = pd.read_csv(assocdb, sep=";")
    evdb = pd.read_csv(evdb, sep=";")
    stadb = pd.read_csv(stadb, sep=";")

    evids = evdb.evid_inst.unique()
    stas = stadb.sta.unique()

    # Loop over available evid_inst
    for evid in evids:
        insacs = inpath_wf + '/sac_' + str(evid) + '/*'
        try:
            st_evid = read(insacs)
        except:
            if verbose: 
                print("event with evid_inst ", evid, "not in inpath_wf")
            continue

        #first loop over stations, to rotate signals
        evtraces = Stream()
        for sta in stas:
            st_sta = st_evid.copy()
            st_sta = st_sta.select(station=sta, component='[E|N]')
            if not st_sta:
                continue
            if len(st_sta) != 2:
                if verbose: 
                    print('evid_inst ' + str(evid) + ' - station ' + sta +
                         ' has wrong n. of chans - will be skipped')
                continue
            #get association info for current evid_inst and sta from assocdb
            assoc_or = assocdb.loc[assocdb['evid_inst'] == evid]
            assoc_or_ev = assoc_or.loc[assoc_or['sta'] == sta]
            baz = assoc_or_ev.iloc[0]['baz']
            st_sta.rotate(method="NE->RT",back_azimuth = baz)
            evtraces.append(st_sta[0])
            evtraces.append(st_sta[1])

        #prepare the output file
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        outf = open(outpath + '/' + str(evid) + '.pkl', 'wb')
        wfpkl = {'signal':{'R':[],'T':[]}, 
                 'noise':{'R':[],'T':[]}, 
                 'stacode':[],
                 'hypodist':[]}
        st_signal_r = Stream()
        st_signal_t = Stream()
        st_noise_r = Stream()
        st_noise_t = Stream()

        #second loop over stations, to process rotated traces
        for tr in evtraces:
            sta = tr.stats['station']
            chan = tr.stats['channel']
            orientation = chan[-1]
            if orientation not in ['R', 'T']:
                continue
            #get association info for current evid and sta from assocdb
            assoc_or = assocdb.loc[assocdb['evid_inst'] == evid]
            assoc_or_ev = assoc_or.loc[assoc_or['sta'] == sta]
            try:
                pickP = assoc_or_ev.loc[assoc_or_ev['phase'] == "P"]
                pickS = assoc_or_ev.loc[assoc_or_ev['phase'] == "S"]
                t0_noise = UTCDateTime(pickP.iloc[0]['picktime']) - tdura -1
                t0_signal = UTCDateTime(pickS.iloc[0]['picktime'])
            except:
                if verbose: 
                    print('evid_inst ' + str(evid) +
                          '- P or S pick missing - skipping station ' + sta)
                continue

            signal = tr.slice(t0_signal, t0_signal + tdura)
            noise = tr.slice(t0_noise, t0_noise + tdura)
            
            for tt in [signal, noise]:
                tt.detrend('demean')
                tt.detrend('linear')

            if orientation == 'R':
                st_signal_r.append(signal)
                st_noise_r.append(noise)
            if orientation == 'T':
                st_signal_t.append(signal)
                st_noise_t.append(noise)
                wfpkl['stacode'].append(sta)
                wfpkl['hypodist'].append(assoc_or_ev.iloc[0]['hypodist'])

        wfpkl['signal']['R'] = st_signal_r
        wfpkl['signal']['T'] = st_signal_t
        wfpkl['noise']['R'] = st_noise_r
        wfpkl['noise']['T'] = st_noise_t
        pickle.dump(wfpkl, outf, protocol=pickle.HIGHEST_PROTOCOL)
        outf.close()

def create_wfpool_mseed(inpath_wf=cfg.WFPATH, outpath=cfg.WFDBPATH, 
                        evdb=cfg.EVDB, stadb=cfg.STADB, assocdb=cfg.ASSOCDB, 
                        tdura=cfg.TDURA, removeresp=True, 
                        respfiles=cfg.RESPPATH, verbose=0):
    
    def clean_trace_metadata(trace):
        trace.stats.processing = []
        trace.stats._format = None
        trace.stats.response = None
        return trace

    assocdb = pd.read_csv(assocdb, sep=";")
    evdb = pd.read_csv(evdb, sep=";")
    stadb = pd.read_csv(stadb, sep=";")

    evids_inst = evdb.evid_inst.unique()
    stas = stadb.sta.unique()

    pltcount = 0

    for evid in evids_inst:
        # skip if the pkl already exists
        outfile = os.path.join(outpath, f"{evid}.pkl")          
        if os.path.isfile(outfile):                             
            if verbose:                                         
                print(f"evid_inst {evid}: pkl file already present, skip.")  
            continue                                           
       

        inwf = os.path.join(inpath_wf, str(evid), '*')
        try:
            st_evid = read(inwf)
        except:
            if verbose:
                print("event with evid_inst ", evid, "not in inpath_wf")
            continue

        if not os.path.exists(outpath):
            os.makedirs(outpath)
        outpath_plt = os.path.join(outpath, 'plots')
        if not os.path.exists(outpath_plt):
            os.makedirs(outpath_plt)

        outf = open(outfile, 'wb')  
        wfpkl = {'signal': {'R': [], 'T': []},
                 'noise': {'R': [], 'T': []},
                 'stacode': [],
                 'hypodist': []}

        st_signal_r, st_signal_t = Stream(), Stream()
        st_noise_r, st_noise_t = Stream(), Stream()

        for sta in stas:
            st_sta_raw = st_evid.select(station=sta, component='[E|N]')
            if not st_sta_raw or len(st_sta_raw) != 2:
                if verbose:
                    print(f"evid_inst {evid} - station {sta} has wrong n. of chans - skipped")
                continue

            # Clean before deepcopy to avoid MemoryError
            for tr in st_sta_raw:
                clean_trace_metadata(tr)
            st_sta = st_sta_raw.copy()

            try:
                assoc_or_ev = assocdb[(assocdb['evid_inst'] == evid) & (assocdb['sta'] == sta)]
                baz = assoc_or_ev.iloc[0]['baz']
            except:
                if verbose:
                    print(f"evid_inst {evid} - station {sta} has no association info - skipped")
                continue

            try:
                pickP = assoc_or_ev[assoc_or_ev['phase'] == "P"]
                pickS = assoc_or_ev[assoc_or_ev['phase'] == "S"]
                t0_noise = UTCDateTime(pickP.iloc[0]['picktime']) - tdura - 1
                t0_signal = UTCDateTime(pickS.iloc[0]['picktime'])
            except:
                if verbose:
                    print(f"evid_inst {evid} - P or S pick missing - skipping station {sta}")
                continue

            signal = st_sta.slice(t0_signal, t0_signal + tdura)
            noise = st_sta.slice(t0_noise, t0_noise + tdura)

            for tt in [signal, noise]:
                tt.detrend('demean')
                tt.detrend('linear')

            if removeresp:
                if not respfiles:
                    if verbose:
                        print(f"no response information for station {sta} - skipped")
                    continue
                inv = read_inventory(os.path.join(respfiles, '*'))
                inv = inv.select(station=sta)
                if len(inv) == 0:
                    if verbose:
                        print(f"response info for station {sta} not found - skipped")
                    continue

                figname = f"{evid}_{st_sta[0].id}_"
                for tt in signal:
                    tt.remove_response(inventory=inv, pre_filt=cfg.PREFILT,
                                       output='VEL', plot=True)
                    plt.savefig(outpath_plt + '/' + figname + str(pltcount)+'.png')
                    plt.close()
                    pltcount+=1
                for tt in noise:
                    tt.remove_response(inventory=inv, pre_filt=cfg.PREFILT,
                                       output='VEL', plot=True)

                    plt.savefig(outpath_plt + '/' + figname + str(pltcount)+'.png')
                    plt.close()
                    pltcount+=1

            try:
                rotated_signal = signal.rotate(method="NE->RT", back_azimuth=baz)
                rotated_noise = noise.rotate(method="NE->RT", back_azimuth=baz)

                if len(rotated_signal) != 2 or len(rotated_noise) != 2:
                    if verbose:
                        print(f"[WARN] Rotation failed for event {evid}, station {sta}.")
                    continue

                signal_t, signal_r = rotated_signal
                noise_t, noise_r = rotated_noise

                st_signal_r.append(signal_r)
                st_noise_r.append(noise_r)
                st_signal_t.append(signal_t)
                st_noise_t.append(noise_t)
                wfpkl['stacode'].append(sta)
                wfpkl['hypodist'].append(assoc_or_ev.iloc[0]['hypodist'])

            except Exception as e:
                if verbose:
                    print(f"[ERROR] Rotation error for event {evid}, station {sta}: {e}")
                continue

        wfpkl['signal']['R'] = st_signal_r
        wfpkl['signal']['T'] = st_signal_t
        wfpkl['noise']['R'] = st_noise_r
        wfpkl['noise']['T'] = st_noise_t
        pickle.dump(wfpkl, outf, protocol=pickle.HIGHEST_PROTOCOL)
        outf.close()

