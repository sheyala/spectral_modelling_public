"""
Constants needed in FAS calculation
"""

import math

PI = math.pi
EE = math.e

RHO = 2.8  # g/cm^3
V_S = 350000  # cm/s

R0 = 100000.  # cm
Rref = 7000000.  # cm

R1 = 4000000.  # cm
R2 = 5000000.  # cm
R3 = 6000000.  # cm
R4 = 8000000.  # cm
R5 = 10000000.  # cm

R_THETA = 0.55
FF = 2.
CSI = 1. / math.sqrt(2.)
LN_MCOST = math.log((R_THETA * FF * CSI) / (4. * PI * RHO * (V_S ** 3) * R0))
MCOST = (R_THETA * FF * CSI) / (4. * PI * RHO * (V_S ** 3) * R0)

QCOST = PI / V_S

AVG_STDROP = 7.3e+06  # g/(cm*s^2)

#for OGS dataset, from Oasis
soil_dict = {'BAD': 'A', 'BOO': 'A', 'BUA': 'A', 'CLUD': 'A', 'FVI': 'A', 
             'GEPF': 'C', 'MLN': 'A', 'MPRI': 'A', 'PLRO': 'A', 'PTCC': 'A', 
             'VINO': 'B', 'ZOU2': 'A', 'ACOM': 'A', 'CIMO': 'A', 
             'COLI': 'A', 'CSM': 'A', 'LSR': 'A', 'PURA': 'C', 
             'STAL': 'B', 'PRED': 'A', 'APGO': 'A', 'CSO': 'A', 'FUSE': 'A'}
