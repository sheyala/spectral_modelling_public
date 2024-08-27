import spectral_modelling.db_create.parsers as parsers
import spectral_modelling.utils.config as cfg
from spectral_modelling.db_create.step1_create_wfpool import create_wfpool
from spectral_modelling.db_create.step2_filter_run import filter_evstas
from spectral_modelling.db_create.step3_save_sptdb import create_sptpkl
import os

def main_parse(sta_filein, staftype, 
               ev_filein, evftype, 
               assoc_filein, assocftype):
    parsers.event_parser(ev_filein, evftype)
    parsers.station_parser(sta_filein, staftype)
    parsers.assoc_parser(assoc_filein, assocftype, stationdb_in=cfg.STADB,
                          eventdb_in=cfg.EVDB)

#TODO fix various cases
def main_steps(step1=True, step2=True, step3=True, format_step1='mseed'):
    if step1:
        create_wfpool(format_step1, verbose=True)
    if step2:
        filter_evstas(verbose=True)
    if step3:
        create_sptpkl()

sta_pathin = cfg.RESPPATH 
ev_pathin = cfg.EVPATH 
assoc_pathin = cfg.EVPATH 

# run code for OGS dataset
if 1:
    main_parse(sta_pathin, "stationxml", ev_pathin, "quakeml", assoc_pathin, "quakeml")
    main_steps(step1=True, step2=True, step3=True)

# # run code for ICELAND dataset
# if 1:
#     main_parse(sta_pathin, "sac", ev_pathin, "csv", assoc_pathin, "csv")
#     main_steps(step1=True, step2=True, step3=True)