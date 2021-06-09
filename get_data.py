import pickle
import scipy.io
import numpy as np
from matplotlib.dates import date2num

"""
A collection of tools to import and standardize data from FVCOM, radar, and ADCP. Assumed formats:
- FVCOM: Pickle (.p) file containing FVCOM data interpolated onto radar grid (regular grid)
- Radar: Pickle (.p) file containing regular grid of data at one level (surface)
- ADCP: MatLab (.mat) file containing time series of velocity components at deployment location
        (top sigma layer, corresponding w sigma = 0.85, is provided)

Provides functions:
    - getFVCOM: grab data from FVCOM file and put it in a dictionary
    - getRadar: grab data from radar file and put it in a dictionary
    - getADCP: grab data from ADCP file and put it in a dictionary

TO-DO:

Lilli Enders (lilli.enders@outlook.com)
June 2021
"""
def getFVCOM(file):
    infile = open(file, 'rb')
    dataIn = pickle.load(infile)
    infile.close()
    data = {}
    keys = ['lon','lat']
    for key in keys:
        data[key] = dataIn[key]
    data['u'] = dataIn['ua']
    data['v'] = dataIn['va']
    data['speed'] = np.sqrt(np.square(data['u']) + np.square(data['v']))
    data['time'] = date2num(dataIn['datetime'])

    startTime = dataIn['datetime'][0]
    endTime = dataIn['datetime'][-1]
    stepTime = dataIn['datetime'][1] - dataIn['datetime'][0]
    data['timeReg'] = date2num(np.arange(startTime, endTime, stepTime))

    return data

def getRadar(file):
    infile = open(file, 'rb')
    dataIn = pickle.load(infile)
    infile.close()
    data = {}
    keys = ['lon','lat']
    for key in keys:
        data[key] = dataIn[key]
    data['u'] = dataIn['ux']
    data['v'] = dataIn['uy']
    data['speed'] = np.sqrt(np.square(data['u'])+np.square(data['v']))
    data['time'] = date2num(dataIn['datetime'])

    startTime = dataIn['datetime'][0]
    endTime = dataIn['datetime'][-1]
    stepTime = dataIn['datetime'][1] - dataIn['datetime'][0]
    data['timeReg'] = date2num(np.arange(startTime, endTime, stepTime))

    return data

def getADCP(file,loc):
    infile = open(file, 'rb')
    dataIn = scipy.io.loadmat(infile)
    infile.close()
    data = {}
    data['matlabtime'] = np.stack(dataIn['matlabtime'])[0]
    data['lon'] = loc[0]
    data['lat'] = loc[1]
    data['u'] = np.stack(dataIn['u'])[0]
    data['v'] = np.stack(dataIn['v'])[0]
    return data

