#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

# import os
# import pandas as pd
# import spectral_modelling.utils.config as cfg


import obspy
from obspy.clients.fdsn.mass_downloader import CircularDomain,\
    Restrictions, MassDownloader

cl_GFZ = obspy.clients.fdsn.Client("GFZ")
starttime=obspy.UTCDateTime("2024-04-23T00:05:00")
endtime=obspy.UTCDateTime("2024-04-23T00:10:00")

#use mass_downloader to get test data for GE.MATE
domain = CircularDomain(latitude=40.65, longitude=16.70,
                        minradius=0.0, maxradius=0.3)  
restrictions = Restrictions(
    starttime=starttime,
    endtime=endtime,
    network="GE",
    station="MATE",
    location="",
    channel="HH*",
    reject_channels_with_gaps=True)
mdl = MassDownloader(providers=[cl_GFZ])
mdl.download(domain, restrictions, mseed_storage='GE_MATE_mseed',
            stationxml_storage='GE_MATE_metadata')
mseed_from_mdownl = obspy.read('GE_MATE_mseed/*')

#use Client to get test data for GE.MATE
mseed_from_client = cl_GFZ.get_waveforms("GE", "MATE", "*", "HH*", 
                                         starttime, endtime)

#compare results from Client and mass_downloader
print(mseed_from_client)
print(mseed_from_mdownl)







# cl_OGS=Client("http://158.110.30.217:8080")
# cl_EIDA=Client("EIDA")
# cl_RER=Client("http://risk.crs.ogs.it:8008")


# #first, load event information from txt file
# eventdb_in=cfg.EVDB
# evdb = pd.read_csv(eventdb_in, sep=";")
# evids = evdb['evid'].tolist()

# #for each event, get event time and request mseeds accordingly
# # for evid in evids:
# evid = evids[0]
# if 0:
#     erow = evdb.loc[evdb['evid'] == evid]
#     etime = erow.iloc[0]['etime']
#     elat = erow.iloc[0]['elat']
#     elon = erow.iloc[0]['elon']
#     evid_inst = erow.iloc[0]['evid_inst']
#     starttime = obspy.UTCDateTime(etime)-60.
#     endtime = obspy.UTCDateTime(etime)+120.
#     print(starttime)
#     print(endtime) 

#     domain = CircularDomain(latitude=44.51, longitude=11.36,
#                             minradius=0.0, maxradius=10.)  

#     restrictions = Restrictions(
#         # Get data for the time window evtime - 60s, evtime + 120s
#         starttime=starttime,
#         endtime=endtime,
#         chunklength_in_sec=180,
#         # Considering the enormous amount of data associated with continuous
#         # requests, you might want to limit the data based on SEED identifiers.
#         # If the location code is specified, the location priority list is not
#         # used; the same is true for the channel argument and priority list.
#         channel_priorities=["HH[ZNE]", "EH[ZNE]", "HN[ZNE]", "HG[ZNE]"],
#         location_priorities=["", "00", "10"],
#         network="OX",
#         station="*",
#         # The typical use case for such a data set are noise correlations where
#         # gaps are dealt with at a later stage.
#         reject_channels_with_gaps=False,
#         # Same is true with the minimum length. All data might be useful.
#         minimum_length=0.0,
#         # Guard against the same station having different names.
#         minimum_interstation_distance_in_m=100.0 )

#     mdl = MassDownloader(providers=[cl_OGS])
#     mdl.download(domain, restrictions, mseed_storage='mdtest',
#                 stationxml_storage='mdtest_xml')




