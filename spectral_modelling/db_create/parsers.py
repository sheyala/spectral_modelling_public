import os
import numpy as np
import pandas as pd
from obspy import read, read_events, read_inventory
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel
from spectral_modelling.utils import config as cfg
from spectral_modelling.utils import utils


# ------------------------------------------------------------------
# utility to parse csv Bulletin file into event dataframe
# ------------------------------------------------------------------
def _read_csv_events(csv_path_in):
    evid = 1
    in_csvs = os.listdir(csv_path_in)
    for csvfile in in_csvs:
        # load each input event csv and keep the first occurrence for each event
        csv_in = pd.read_csv(os.path.join(csv_path_in, csvfile), sep=";")
        csv_in = csv_in.drop_duplicates("ID")
        csv_in = csv_in.sort_values("Date")
        csv_in = csv_in.reset_index()
        df_out = pd.DataFrame(columns=["evid_inst", "etime", "elat",
                                       "elon", "edepth", "ml", "mw"])
        for i in range(len(csv_in)):
            row = csv_in.loc[i]
            mw = "{:.2f}".format(utils.ml_to_mag_munafo(float(row["Mag"])))
            df_out.loc[evid] = pd.Series({
                "evid_inst": row["ID"],
                "etime":     row["Date"] + "T" + row["Time"] + "Z",
                "elat":      float(row["Lat"][:-1]),
                "elon":      float(row["Lon"][:-1]),
                "edepth":    "{:.3f}".format(row["Depth"]),
                "ml":        row["Mag"],
                "mw":        mw
            })
            evid += 1
    df_out.index.name = "evid"
    return df_out


# ------------------------------------------------------------------
# utility to parse quakeml file into event dataframe
# ------------------------------------------------------------------
def _read_qml_events(qml_path_in):
    from obspy import read_events

    df_out = pd.DataFrame(columns=["evid_inst", "etime", "elat",
                                   "elon", "edepth", "ml", "mw"])
    evid = 1
    in_qmls = os.listdir(cfg.EVPATH)
    for qmlfile in in_qmls:
        qml_in = read_events(os.path.join(qml_path_in, qmlfile))
        for ev in qml_in.events:
            # institution evid (remove prefix if present)
            evid_inst = ev.resource_id.id
            try:
                evid_inst = int(evid_inst.split("?")[1].split("=")[1])
            except Exception:
                evid_inst = int(str(ev.resource_id).split('/')[-1])

            origin   = ev.preferred_origin()
            edepth   = origin.depth / 1000.0
            elat     = origin.latitude
            elon     = origin.longitude
            etime    = origin.time

            magnitude = ev.preferred_magnitude()
            magtype   = magnitude.magnitude_type

            if magtype.lower() == "md":
                print(f"[IGNORATO] Evento {evid_inst} con magnitudo tipo '{magtype}' ignorato.")
                continue

            if magtype.lower() == "ml":
                ml = magnitude.mag
                mw = utils.ml_to_mag_munafo(ml)
            elif magtype.lower() == "mw":
                mw = magnitude.mag
                ml = utils.mag_to_ml_munafo(mw)
            else:
                print(f"[IGNORATO] Evento {evid_inst} con magnitudo tipo non gestito '{magtype}' ignorato.")
                continue

            df_out.loc[evid] = pd.Series({
                "evid_inst": evid_inst,
                "etime":     etime,
                "elat":      elat,
                "elon":      elon,
                "edepth":    edepth,
                "ml":        f"{ml:.2f}",
                "mw":        f"{mw:.2f}"
            })
            evid += 1

    df_out.index.name = "evid"
    return df_out


# ------------------------------------------------------------------
# utility to parse sac files into station dataframe
# ------------------------------------------------------------------
def _read_sac_stations(sac_path_in):
    df_out = pd.DataFrame(columns=["net", "sta", "slat", "slon", "selev"])
    in_sacs  = os.listdir(sac_path_in)
    oldstas  = []
    i = 1
    for sacfolder in in_sacs:
        st_sac = read(os.path.join(sac_path_in, sacfolder, "*E.sac"))
        for tr in st_sac:
            sacinfo = tr.stats.sac
            if sacinfo["kstnm"] not in oldstas:
                df_out.loc[i] = pd.Series({
                    "net":   cfg.NETCODE,
                    "sta":   sacinfo["kstnm"],
                    "slat":  float(sacinfo["stla"]),
                    "slon":  float(sacinfo["stlo"]),
                    "selev": float(sacinfo["stel"])
                })
                oldstas.append(sacinfo["kstnm"])
                i += 1
    df_out = df_out.sort_values("sta")
    return df_out


# ------------------------------------------------------------------
# utility to parse stationxml files into station dataframe
# ------------------------------------------------------------------
def _read_xml_stations(xml_path_in):
    df_out = pd.DataFrame(columns=["net", "sta", "slat", "slon", "selev"])
    in_xmls = os.listdir(xml_path_in)
    oldstas = []
    i = 1
    for in_xml in in_xmls:
        st_xml = read_inventory(os.path.join(xml_path_in, in_xml))
        for net in st_xml.networks:
            for sta in net.stations:
                if sta.code not in oldstas:
                    df_out.loc[i] = pd.Series({
                        "net":   net.code,
                        "sta":   sta.code,
                        "slat":  sta.latitude,
                        "slon":  sta.longitude,
                        "selev": sta.elevation
                    })
                    oldstas.append(sta.code)
                    i += 1
    df_out = df_out.sort_values("sta")
    return df_out


# ------------------------------------------------------------------
# utility to parse csv Bulletin file into associations dataframe
# ------------------------------------------------------------------
def _read_csv_assoc(csv_path_in, station_db_in, event_db_in):
    stadb = pd.read_csv(station_db_in, sep=";")
    evdb  = pd.read_csv(event_db_in, sep=";")
    evids       = evdb["evid"].tolist()
    evids_inst  = evdb["evid_inst"].tolist()

    model = TauPyModel(model="iasp91")

    aid = 1
    in_csvs = os.listdir(csv_path_in)
    df_out = pd.DataFrame(columns=["evid", "evid_inst", "sta", "phase",
                                   "picktime", "baz", "takeoff", "hypodist"])
    for csvfile in in_csvs:
        csv_in = pd.read_csv(os.path.join(csv_path_in, csvfile), sep=";")
        csv_in = csv_in.sort_values("Date")

        orids = csv_in.ID.unique()
        for orid in orids:
            if orid not in evids_inst:
                continue
            evid = evids[evids_inst.index(orid)]

            csv_orid = csv_in.loc[csv_in["ID"] == orid]
            stas = np.sort(csv_orid.Station.unique())

            for sta in stas:
                csv_sta = csv_orid.loc[csv_orid["Station"] == sta]
                csv_sta = csv_sta.sort_values("Time_pick").reset_index()

                stameta  = stadb.loc[stadb["sta"] == sta].iloc[0]
                slat, slon = stameta["slat"], stameta["slon"]
                elat, elon = float(csv_sta.iloc[0]["Lat"][:-1]), float(csv_sta.iloc[0]["Lon"][:-1])
                edepth     = csv_sta.iloc[0]["Depth"]

                epidist, baz, _az = gps2dist_azimuth(slat, slon, elat, elon)
                baz       = "{:.2f}".format(baz)
                hypodist  = "{:.3f}".format(np.sqrt(edepth**2 + (epidist / 1000.0) ** 2))

                dist_deg  = locations2degrees(slat, slon, elat, elon)
                arrivals  = model.get_travel_times(source_depth_in_km=edepth,
                                                   distance_in_degree=dist_deg)
                takeoff   = "{:.2f}".format(arrivals[0].takeoff_angle)

                for _, row in csv_sta.iterrows():
                    df_out.loc[aid] = pd.Series({
                        "evid":      evid,
                        "evid_inst": row["ID"],
                        "sta":       row["Station"],
                        "phase":     row["Phase"],
                        "picktime":  row["Time_pick"],
                        "baz":       baz,
                        "takeoff":   takeoff,
                        "hypodist":  hypodist
                    })
                    aid += 1
    return df_out


# ------------------------------------------------------------------
# utility to parse quakeml file into associations dataframe
# ------------------------------------------------------------------
def _read_qml_assoc(xml_path_in, station_db_in, event_db_in):
    stadb       = pd.read_csv(station_db_in, sep=";")
    stas        = stadb["sta"].tolist()
    evdb        = pd.read_csv(event_db_in, sep=";")
    evids       = evdb["evid"].tolist()
    evids_inst  = evdb["evid_inst"].tolist()

    model = TauPyModel(model="iasp91")

    aid = 1
    df_out = pd.DataFrame(columns=["evid", "evid_inst", "sta", "phase",
                                   "picktime", "baz", "takeoff", "hypodist"])
    for xmlfile in os.listdir(xml_path_in):
        qml_in = read_events(os.path.join(xml_path_in, xmlfile))
        for ev in qml_in.events:
            evid_inst = ev.resource_id.id
            try:
                evid_inst = int(evid_inst.split("?")[1].split("=")[1])
            except Exception:
                evid_inst = int(str(ev.resource_id).split('/')[-1])

            if evid_inst not in evids_inst:
                continue
            evid = evids[evids_inst.index(evid_inst)]

            picks = ev.picks
            for sta in stas:
                picks_sta_p = [p for p in picks if p.waveform_id.station_code == sta and p.phase_hint == "P"]
                picks_sta_s = [p for p in picks if p.waveform_id.station_code == sta and p.phase_hint == "S"]
                if not picks_sta_p and not picks_sta_s:
                    continue

                stameta          = stadb.loc[stadb["sta"] == sta].iloc[0]
                slat, slon       = stameta["slat"], stameta["slon"]
                origin           = ev.preferred_origin()
                elat, elon, edepth = origin.latitude, origin.longitude, origin.depth / 1000.0

                epidist, baz, _az = gps2dist_azimuth(slat, slon, elat, elon)
                baz       = "{:.2f}".format(baz)
                hypodist  = "{:.3f}".format(np.sqrt(edepth**2 + (epidist / 1000.0) ** 2))

                dist_deg   = locations2degrees(slat, slon, elat, elon)
                arrivals   = model.get_travel_times(source_depth_in_km=edepth,
                                                    distance_in_degree=dist_deg)
                takeoff    = "{:.2f}".format(arrivals[0].takeoff_angle)

                # pick P
                if picks_sta_p:
                    if len(picks_sta_p) > 1:
                        picks_sta_p = [p for p in picks_sta_p if p.evaluation_mode == "manual"]
                    if len(picks_sta_p) > 1:
                        uncert = [p.time_errors.uncertainty for p in picks_sta_p if p.time_errors]
                        if uncert:
                            picks_sta_p = [picks_sta_p[uncert.index(min(uncert))]]
                    if len(picks_sta_p) == 1:
                        df_out.loc[aid] = pd.Series({
                            "evid":      evid,
                            "evid_inst": evid_inst,
                            "sta":       sta,
                            "phase":     "P",
                            "picktime":  picks_sta_p[0].time,
                            "baz":       baz,
                            "takeoff":   takeoff,
                            "hypodist":  hypodist
                        })
                        aid += 1

                # pick S
                if picks_sta_s:
                    if len(picks_sta_s) > 1:
                        picks_sta_s = [p for p in picks_sta_s if p.evaluation_mode == "manual"]
                    if len(picks_sta_s) > 1:
                        uncert = [p.time_errors.uncertainty for p in picks_sta_s if p.time_errors]
                        if uncert:
                            picks_sta_s = [picks_sta_s[uncert.index(min(uncert))]]
                    if len(picks_sta_s) == 1:
                        df_out.loc[aid] = pd.Series({
                            "evid":      evid,
                            "evid_inst": evid_inst,
                            "sta":       sta,
                            "phase":     "S",
                            "picktime":  picks_sta_s[0].time,
                            "baz":       baz,
                            "takeoff":   takeoff,
                            "hypodist":  hypodist
                        })
                        aid += 1
    return df_out


# ------------------------------------------------------------------
# utility to create csv station list from input
# ------------------------------------------------------------------
def station_parser(filein, filetype):
    # -------- SKIP if the file already exists --------
    if os.path.isfile(cfg.STADB):
        print(f"{cfg.STADB} già presente – salta station_parser.")
        return
    # --------------------------------------------
    match filetype:
        case "csv":
            pass  
        case "stationxml":
            dfout = _read_xml_stations(filein)
        case "sac":
            dfout = _read_sac_stations(filein)
        case _:
            print("Filetype not recognized")
            return
    dfout.to_csv(cfg.STADB, sep=";", index=False)


# ------------------------------------------------------------------
# utility to create csv event list from input
# ------------------------------------------------------------------
def event_parser(filein, filetype):
    # -------- SKIP if the file already exists --------
    if os.path.isfile(cfg.EVDB):
        print(f"{cfg.EVDB} già presente – salta event_parser.")
        return
    # --------------------------------------------
    match filetype:
        case "csv":
            dfout = _read_csv_events(filein)
        case "quakeml":
            dfout = _read_qml_events(filein)
        case _:
            print("Filetype not recognized")
            return
    dfout.to_csv(cfg.EVDB, sep=";")


# ------------------------------------------------------------------
# utility to create csv associations list from input
# ------------------------------------------------------------------
def assoc_parser(filein, filetype, stationdb_in=None, eventdb_in=None):
    # -------- SKIP if the file already exists --------
    if os.path.isfile(cfg.ASSOCDB):
        print(f"{cfg.ASSOCDB} già presente – salta assoc_parser.")
        return
    # --------------------------------------------
    match filetype:
        case "csv":
            if not stationdb_in or not eventdb_in:
                print("stationdb_in ed eventdb_in devono essere forniti per il formato csv.")
                return
            dfout = _read_csv_assoc(filein, stationdb_in, eventdb_in)
        case "quakeml":
            if not stationdb_in or not eventdb_in:
                print("stationdb_in ed eventdb_in devono essere forniti per il formato quakeml.")
                return
            dfout = _read_qml_assoc(filein, stationdb_in, eventdb_in)
        case _:
            print("Filetype not recognized")
            return
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
