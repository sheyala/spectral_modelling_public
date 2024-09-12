# -----
# Configuration File for db_create and inversion
# -----
# BASEPATH = home folder for calculations
# DATAPATH = folder for data, metadata and databases
# WFPATH   = folder for original waveform files (input for step 1)
# WFDBPATH = folder for wafeform databases (input for step 2 and step 3)
# SPTDBPATH= folder containing pickle spectrum databases
#            (input for inversion)
# DBNAME   = name of the "database version", e.g. NEI for the 
#            North-East Italy dataset
# RUNNAME  = name of the run, uniquely assigned to the set of data constraints
#            used in Step 2

import numpy as np
import os

BASEPATH = "/home/laura/Desktop/per_shakelab/spectral_modelling_public/spectral_modelling"
DATAPATH = os.path.join(BASEPATH, "data") 

#OGS dataset case
DBNAME = "OGS2019"
DBPATH = os.path.join(DATAPATH, DBNAME)
WFPATH = os.path.join(DBPATH, "mseeds")
RESPPATH = os.path.join(DBPATH, "stas_xml")
EVPATH = os.path.join(DBPATH, "event_qml")

# # ICELAND dataset case
# DBNAME = "REYK"
# DBPATH = os.path.join(DATAPATH, DBNAME)
# WFPATH = os.path.join(DBPATH, "event_sac")
# RESPPATH = os.path.join(DBPATH, "event_sac")
# EVPATH = os.path.join(DBPATH, "assoc_bul")
# NETCODE = '7E'


WFDBPATH = os.path.join(DBPATH, "wfdb")
SPTDBPATH = os.path.join(DBPATH, "sptdb")
RUNNAME = "run1"
RUNPATH = os.path.join(DBPATH, RUNNAME + "_output")

# WFPATH = os.path.join(DATAPATH, "event_sac")
# DBNAME = "REYK"
# DBNAME = "INGV2"

#-----------------------------------------------
#configuration for PARSER and STEP 1 
#parse the SAC db and save cut waveform db for each event
#-----------------------------------------------
EVDB = os.path.join(DBPATH, DBNAME + "_events_all.txt")
STADB = os.path.join(DBPATH, DBNAME + "_stas_all.txt")
ASSOCDB = os.path.join(DBPATH, DBNAME + "_assoc_all.txt")
TDURA = 30 # duration of the signal (from S pick) 
           # and noise (pre P pick) windows
PREFILT = [0.05, 0.1, 95., 100.] #optional pre-filter for response removal

#-----------------------------------------------
#configuration for STEP 2 
#apply constraints on the dataset and filter recursively
#-----------------------------------------------
SNRMIN = 3.
F0_SNR = 2.
F1_SNR = 15.
FMINSTOP = 10
FMAXSTOP = 10
MIN_EVS = 3
MIN_STAS = 3
MIN_HYPODIST = 5.
MAX_HYPODIST = 100.

#-----------------------------------------------
#configuration for STEP 3 - 
#save spectra databases (one for each of R/T components, signal/noise)
#-----------------------------------------------
EVDB_FOR_INVERSION = os.path.join(DBPATH, DBNAME + "_events_" + RUNNAME +
                                  ".txt")
STADB_FOR_INVERSION = os.path.join(DBPATH, DBNAME + "_stas_" + RUNNAME +
                                  ".txt")
FREQLIMS = [0.5, 25.]

#-----------------------------------------------
#configuration for INVERSION (and for utils)
#-----------------------------------------------
#tbd
COMPNAME = dict(S="_signal_", 
                N="_noise_"
                )

FREQLIST = np.logspace(np.log10(FREQLIMS[0]), np.log10(FREQLIMS[1]), 
                       num=30, endpoint=True)
FREQS = np.array(FREQLIST, dtype=np.double)
KO_B = 40.
SYNTHPARS_PATH = os.path.join(DBPATH, "synthpars")
RUNPLOT_PATH = os.path.join(RUNPATH, "plots")

# LOG_FOLDER = os.path.join(DATAPATH, DBNAME + "_logs")

# lists of recognized minimazion methods
MINIMIZE_METHODS_BOUNDS = ['Powell', 
                           'L-BFGS-B', 
                           'TNC', 
                           'trust-constr',
                           'SLSQP']

MINIMIZE_METHODS_NOBOUNDS = ['Nelder-Mead', 
                             'Powell', 
                             'CG', 
                             'BFGS',
                             'Newton-CG', 
                             'L-BFGS-B', 
                             'TNC', 
                             'COBYLA',
                             'SLSQP', 
                             'trust-constr', 
                             'dogleg', 
                             'trust-ncg',
                             'trust-exact', 
                             'trust-krylov']