import os
import pandas as pd
from spectral_modelling.utils import config as cfg
from spectral_modelling.utils import utils
from obspy import read, read_events, read_inventory
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees
import numpy as np
from obspy.taup import TauPyModel


#utility to parse csv Bulletin file into event dataframe
def _read_csv_events(csv_path_in):
    evid = 1
    in_csvs = os.listdir(csv_path_in)
    for csvfile in in_csvs:
        #load each input event csv and keep the first occurence for each event
        csv_in = pd.read_csv(os.path.join(csv_path_in, csvfile), sep=";")
        csv_in = csv_in.drop_duplicates("ID")
        csv_in = csv_in.sort_values("Date")
        csv_in = csv_in.reset_index()
        df_out = pd.DataFrame(columns=["evid_inst","etime","elat",
                                    "elon", "edepth", "ml", "mw"])
        for i in range(len(csv_in)):
            row = csv_in.loc[i]
            mw = str("{:.2f}".format(utils.ml_to_mag_munafo(float(row["Mag"]))))
            df_out.loc[evid] = pd.Series({"evid_inst":row["ID"],
                                        "etime":row["Date"]+"T"+row["Time"]+"Z",
                                        "elat":float(row["Lat"][:-1]),
                                        "elon":float(row["Lon"][:-1]),
                                        "edepth":str("{:.3f}".format(row["Depth"])),
                                        "ml":row["Mag"],
                                        "mw":mw})
            evid +=1
    df_out.index.name = "evid"
    return df_out


#utility to parse quakeml file into event dataframe
# NOTE the institution evid is kept for reference as evid_inst
def _read_qml_events(qml_path_in):
    #prepare output dataframe
    df_out = pd.DataFrame(columns=["evid_inst","etime","elat",
                                   "elon", "edepth", "ml", "mw"])

    evid = 1
    #load and process each input event quakeml 
    in_qmls = os.listdir(qml_path_in)
    for qmlfile in in_qmls:
        qml_in = read_events(os.path.join(qml_path_in, qmlfile))
        events = qml_in.events
        for ev in events:
            #get institution evid, remove prefix
            evid_inst = ev.resource_id.id
            try:
                evid_inst = int(evid_inst.split("?")[1].split("=")[1])
            except:
                evid_inst = int(str(ev.resource_id).split('/')[-1])
            #get origin information
            origin = ev.preferred_origin()
            edepth = origin.depth/1000. #convert to CGS
            elat = origin.latitude
            elon = origin.longitude
            etime = origin.time

            #get magniutde information
            magnitude = ev.preferred_magnitude()
            magtype = magnitude.magnitude_type
            if magtype not in ['ML', 'Ml', 'ml', 'MW', 'Mw', 'mw']:
                raise TypeError('Magnitude type must be either MW or ML, not %s' % magtype)
                continue
            if magtype in['ML', 'Ml', 'ml']:
                ml = magnitude.mag
                mw = utils.ml_to_mag_munafo(ml)
            if magtype in ['MW', 'Mw', 'mw']:
                mw = magnitude.mag
                ml = utils.mag_to_ml_munafo(mw)

            df_out.loc[evid] = pd.Series({"evid_inst":evid_inst,
                                        "etime":etime,
                                        "elat":elat,
                                        "elon":elon,
                                        "edepth":edepth,
                                        "ml":str("{:.2f}".format(ml)),
                                        "mw":str("{:.2f}".format(mw))})
            evid +=1
    df_out.index.name = "evid"
    return df_out




#utility to parse sac files into station dataframe
def _read_sac_stations(sac_path_in):
    df_out = pd.DataFrame(columns=["net","sta","slat","slon", "selev"])
    in_sacs = os.listdir(sac_path_in)
    oldstas = []
    i = 1
    for in_sac in in_sacs:
        sacpath = os.path.join(sac_path_in, in_sac)
        st_sac = read(sacpath + "/*E.sac")
        for tr in st_sac:
            sacinfo = tr.stats.sac
            #only keep each station once
            if sacinfo["kstnm"] not in oldstas:
                df_out.loc[i] = pd.Series({"net":cfg.NETCODE,
                                           "sta":sacinfo["kstnm"],
                                           "slat":float(str(sacinfo["stla"])),
                                           "slon":float(str(sacinfo["stlo"])),
                                           "selev":float(str(sacinfo["stel"]))})
                i += 1
                oldstas.append(sacinfo["kstnm"])
    df_out = df_out.sort_values("sta")
    return df_out


#utility to parse stationxml files into station dataframe
def _read_xml_stations(xml_path_in):
    df_out = pd.DataFrame(columns=["net","sta","slat","slon", "selev"])
    in_xmls = os.listdir(xml_path_in)
    oldstas = []
    i = 1
    for in_xml in in_xmls:
        st_xml = read_inventory(os.path.join(xml_path_in, in_xml))
        networks = st_xml.networks
        for net in networks:
            stations = net.stations
            for sta in stations:
                if sta.code not in oldstas:
                    df_out.loc[i] = pd.Series({"net":net.code,
                                               "sta":sta.code,
                                               "slat":sta.latitude,
                                               "slon":sta.longitude,
                                               "selev":sta.elevation})
                    i += 1
                    oldstas.append(sta.code)

    df_out = df_out.sort_values("sta")
    return df_out




#utility to parse csv Bulletin file into associations dataframe
def _read_csv_assoc(csv_path_in, station_db_in, event_db_in):
    #load station info
    stadb = pd.read_csv(station_db_in, sep=";")
    #load event info
    evdb = pd.read_csv(event_db_in, sep=";")
    evids = evdb['evid'].tolist()
    evids_inst = evdb['evid_inst'].tolist()

    model = TauPyModel(model="iasp91")

    #load all association csvs
    aid = 1
    in_csvs = os.listdir(csv_path_in)
    for csvfile in in_csvs:
        csv_in = pd.read_csv(os.path.join(csv_path_in, csvfile), sep=";")
        csv_in = csv_in.sort_values("Date")
        df_out = pd.DataFrame(columns=["evid","evid_inst","sta","phase",
                                       "picktime","baz","takeoff","hypodist"])
        orids = csv_in.ID.unique()
        for orid in orids:
            if not orid in evids_inst:
                continue
            evid = evids[evids_inst.index(orid)]
            #select info for current orid and get all associated stations
            csv_orid = csv_in.loc[csv_in["ID"] == orid]
            stas = csv_orid.Station.unique()
            stas.sort()
            for sta in stas:
                #select info for current station
                csv_sta = csv_orid.loc[csv_orid["Station"] == sta]
                csv_sta = csv_sta.sort_values("Time_pick")
                csv_sta = csv_sta.reset_index()

                #calculate distance and baz (same for both picks)
                stameta = stadb.loc[stadb['sta'] == sta].iloc[0]
                slat = stameta['slat']
                slon = stameta['slon']
                elat = float(csv_sta.iloc[0]["Lat"][:-1])
                elon = float(csv_sta.iloc[0]["Lon"][:-1])
                edepth = csv_sta.iloc[0]["Depth"]
                epidist, baz, az = gps2dist_azimuth(slat, slon, elat, elon)
                baz = str("{:.2f}".format(baz))
                hypodist = np.sqrt(edepth**2 + (epidist/1000.)**2)
                hypodist = str("{:.3f}".format(hypodist))

                dist_in_deg = locations2degrees(slat, slon, elat, elon)
                arrivals = model.get_travel_times(source_depth_in_km=edepth, 
                                              distance_in_degree=dist_in_deg)
                takeoff_angle = str("{:.2f}".format(arrivals[0].takeoff_angle))

                #store info in df_out
                for i in range(len(csv_sta)):
                    row = csv_sta.iloc[i]
                    elat = float(row["Lat"][:-1])
                    elon = float(row["Lon"][:-1])
                    df_out.loc[aid] = pd.Series({"evid":evid,
                                                "evid_inst":row["ID"],
                                                "sta":row["Station"],
                                                "phase":row["Phase"],
                                                "picktime":row["Time_pick"],
                                                "baz":baz,
                                                "takeoff":takeoff_angle,
                                                "hypodist":hypodist})
                    aid +=1
    return df_out


#utility to parse quakeml file into associations dataframe
def _read_qml_assoc(xml_path_in, station_db_in, event_db_in):
    df_out = pd.DataFrame(columns=["evid","evid_inst","sta","phase","picktime",
                                   "baz","takeoff","hypodist"])
    #load station info
    stadb = pd.read_csv(station_db_in, sep=";")
    stas = stadb['sta'].tolist()

    #load event info
    evdb = pd.read_csv(event_db_in, sep=";")
    evids = evdb['evid'].tolist()
    evids_inst = evdb['evid_inst'].tolist()

    model = TauPyModel(model="iasp91")

    #load input event quakemls 
    aid = 1
    in_xmls = os.listdir(xml_path_in)
    for xmlfile in in_xmls:
        qml_in = read_events(os.path.join(xml_path_in, xmlfile))
        events = qml_in.events

        for ev in events:
            #get institution evid, remove prefix
            evid_inst = ev.resource_id.id
            try:
                evid_inst = int(evid_inst.split("?")[1].split("=")[1])
            except:
                evid_inst = int(str(ev.resource_id).split('/')[-1])
            if evid_inst not in evids_inst:
                continue
            evid = evids[evids_inst.index(evid_inst)]

            picks = ev.picks
            for sta in stas:
                #select info for current station
                picks_sta_p = [p for p in picks if \
                            (p.waveform_id.station_code == sta) and \
                            (p.phase_hint == 'P')]
                picks_sta_s = [p for p in picks if \
                            (p.waveform_id.station_code == sta) and \
                            (p.phase_hint == 'S')]
                if not picks_sta_s and not picks_sta_p:
                    continue
                
                #calculate distance, baz and takeoff (same for both picks)
                stameta = stadb.loc[stadb['sta'] == sta].iloc[0]
                slat = stameta['slat']
                slon = stameta['slon']
                origin = ev.preferred_origin()
                elat = origin.latitude
                elon = origin.longitude
                edepth = origin.depth/1000. #convert to km
                epidist, baz, az = gps2dist_azimuth(slat, slon, elat, elon)
                baz = str("{:.2f}".format(baz))
                hypodist = np.sqrt(edepth**2 + (epidist/1000.)**2) # in km
                hypodist = str("{:.3f}".format(hypodist))

                dist_in_deg = locations2degrees(slat, slon, elat, elon)
                arrivals = model.get_travel_times(source_depth_in_km=edepth, 
                                              distance_in_degree=dist_in_deg)
                takeoff_angle = str("{:.2f}".format(arrivals[0].takeoff_angle)) 


                #extract best P pick
                if picks_sta_p:
                    #if multiple picks are present, only keep manual ones
                    if len(picks_sta_p) > 1:
                        picks_sta_p = [p for p in picks_sta_p if p.evaluation_mode \
                                    == 'manual']
                    #if multiple picks are present, keep the one with minimum 
                    #associated time uncertainty
                    if len(picks_sta_p) > 1:
                        timeuncert = [p.time_errors.uncertainty for p in \
                                    picks_sta_p]
                        picks_sta_p = \
                            [picks_sta_p[timeuncert.index(min(timeuncert))]]
                    #if multiple picks are still present, raise error
                    if len(picks_sta_p) > 1:   
                        raise TypeError('Cannot disambiguate '+ 
                                        str(len(picks_sta_p)) +
                                        ' P picks for event ' + 
                                        str(evid_inst) + 
                                        ' and station ' + sta)
                    else:
                        #store info in df_out
                        df_out.loc[aid] = pd.Series({"evid":evid,
                                                    "evid_inst":evid_inst,
                                                    "sta":sta,
                                                    "phase":'P',
                                                    "picktime":picks_sta_p[0].time,
                                                    "baz":baz,
                                                    "takeoff":takeoff_angle,
                                                    "hypodist":hypodist})
                        aid +=1

                #extract best S pick
                if picks_sta_s:
                    #if multiple picks are present, only keep manual ones
                    if len(picks_sta_s) > 1:
                        picks_sta_s = [p for p in picks_sta_s if p.evaluation_mode \
                                    == 'manual']
                    #if multiple picks are present, keep the one with minimum 
                    #associated time uncertainty
                    if len(picks_sta_s) > 1:
                        timeuncert = [p.time_errors.uncertainty for p in \
                                    picks_sta_s]
                        picks_sta_s = \
                            [picks_sta_s[timeuncert.index(min(timeuncert))]]
                    #if multiple picks are still present, raise error
                    if len(picks_sta_s) > 1:   
                        raise TypeError('Cannot disambiguate '+ 
                                        str(len(picks_sta_s)) +
                                        ' S picks for event ' + 
                                        str(evid_inst) + 
                                        ' and station ' + sta)
                    else:
                        #store info in df_out
                        df_out.loc[aid] = pd.Series({"evid":evid,
                                                    "evid_inst":evid_inst,
                                                    "sta":sta,
                                                    "phase":'S',
                                                    "picktime":picks_sta_s[0].time,
                                                    "baz":baz,
                                                    "takeoff":takeoff_angle,
                                                    "hypodist":hypodist})
                        aid +=1

    return df_out



#utility to create csv station list from input
def station_parser(filein, filetype):
    match filetype:
        case "csv":
            pass
        case "stationxml":
            dfout =_read_xml_stations(filein)
        case "sac":
            dfout =_read_sac_stations(filein)
        case _: 
            print("Filetype not recognized")
            return None
    dfout.to_csv(cfg.STADB, sep=";", index=False)


#utility to create csv event list from input
def event_parser(filein, filetype):
    match filetype:
        case "csv":
            dfout = _read_csv_events(filein)
        case "quakeml":
            dfout = _read_qml_events(filein)
        case _: 
            print("Filetype not recognized")
            return None
    dfout.to_csv(cfg.EVDB, sep=";")


#utility to create csv associations list from input
def assoc_parser(filein, filetype, stationdb_in=None, eventdb_in=None):
    match filetype:
        case "csv":
            if not stationdb_in:
                print("stationdb_in must be provided for the csv Bulletin case")
                return None 
            if not eventdb_in:
                print("event_db_in must be provided for the csv Bulletin case")
                return None 
            dfout = _read_csv_assoc(filein, stationdb_in, eventdb_in)
        case "quakeml":
            if not stationdb_in:
                print("stationdb_in must be provided for the quakeml case")
                return None 
            if not eventdb_in:
                print("event_db_in must be provided for the quakeml case")
                return None 
            dfout = _read_qml_assoc(filein, stationdb_in, eventdb_in)
        case _: 
            print("Filetype not recognized")
            return None
    dfout.to_csv(cfg.ASSOCDB, sep=";", index=False)


# TODO wrap event parser so that the existing event list is updated if new
# TODO orid is available, instead of created anew
def create_eventdb():
    pass


# TODO wrap station parser so that the existing station list is updated if new
# TODO sta is available, instead of created anew
def create_stationdb():
    pass


# TODO wrap assoc parser so that the existing assoc list is updated if new
# TODO association is available, instead of created anew
def create_assocdb():
    pass