import pickle
import datetime
import numpy as np
import utide as ut
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import pandas as pd
import csv
from datetime import timedelta

# 2D FVCOM data from Jeremy
f = 'Data/FVCOM_2D_Data_2020.p'
infile = open(f, 'rb')
FVCOM = pickle.load(infile,encoding='latin1')
infile.close()

# Name of output file
svname = 'MP_Stations_UTide_June21'
date_cols = ["deploy_date_time","recover_date_time"]
MPdf = pd.read_csv ('MP_Stations.csv',parse_dates=date_cols)

# Get datetimes for u-tide
startTime = datetime.datetime(2020, 1, 1,0,0,0,0)
nextTime = datetime.datetime(2020,1,1,0,10,0)
stepTime =  nextTime - startTime
endTime = datetime.datetime(2021, 1, 1,0,0,0,0)
time = date2num(np.arange(startTime, endTime, stepTime))


# Just organizing ourselves a lil bit
dataDict = FVCOM['data']
outdf = pd.DataFrame({'deploy_year' : [],'station' : [],'mean_lat' : [],
                     'mean_lon' : [],'time' : [],'el' : [],'sl' : []})
# Loop through all keys, get lat and lon at key location, perform utide to get constituents

for i in range(len(MPdf)):
    # Find index of station name
    station = MPdf['station'][i]
    deploy_year = MPdf['deploy_year'][i]
    lat = MPdf['mean_lat'][i]
    lon = MPdf['mean_long'][i]
    outTime = np.arange(MPdf['deploy_date_time'][i],MPdf['recover_date_time'][i],timedelta(hours=1))
    outTimeNum = date2num(outTime)
    el = dataDict[MPdf['station'][i]]['el']
    mwl = dataDict[MPdf['station'][i]]['MWL']
    coefTemp = ut.solve(time, el, lat=dataDict[MPdf['station'][i]]['lat'], conf_int='linear', epoch='python')
    ut_slTemp = ut.reconstruct(outTimeNum, coefTemp, epoch='python')
    for j in range(len(outTime)):
        dftemp = pd.DataFrame(data={'deploy_year' : [deploy_year],'station' : [station],'mean_lat' : [lat],
                     'mean_lon' : [lon],'time' : [outTime[j]],'el' : [ut_slTemp['h'][j]],'sl' : [ut_slTemp['h'][j]+ mwl]})
        outdf = outdf.append(dftemp)
    print(station + 'done')

compression_opts = dict(method='zip',
                        archive_name='MPStations.csv')
outdf.to_csv('MPStations.zip', index=False, compression=compression_opts)