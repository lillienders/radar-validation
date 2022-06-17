from val_plots import *
import numpy as np


"""
Calculations used for radar validation

Provides functions:
    - vorticity: takes in a dictionary (from getFVCOM() or getRadar() functions in get_data), calcuate 
    - divergence: grab u-tide coefficients and use them to reconstruct time-series at each location
    - sshgrad: calculates spatial gradient of sea surface height (SSH) between grid points
    - residuals: from u_tide results, get residuals and calculate vorticity/divergence
Lilli Enders (lilli.enders@outlook.com)
June 2022
"""
def vorticity(dataDict,timeidx=0,plot = False):
    vort = np.empty((len(dataDict['time']),len(dataDict['lat']), len(dataDict['lon'])))
    vort[:] = np.nan
    for t in range (1,len(dataDict['time'])):
        for i in range(1,len(dataDict['lon'])-1):
            for j in range(1,len(dataDict['lat'])-1):
                dvdx = (dataDict['u'][t, j + 1, i + 1] - dataDict['u'][t, j - 1, i - 1]) / 2 * (dataDict['lon'][1] - dataDict['lon'][0])
                dudy = (dataDict['v'][t, j + 1, i + 1] - dataDict['v'][t, j - 1, i - 1]) / 2 * (dataDict['lat'][1] - dataDict['lat'][0])
                vort[t,j, i] = dvdx - dudy
    dataDict['vort'] = vort
    if plot == True:
        bb = [dataDict['lon'][0], dataDict['lon'][-1], dataDict['lat'][0], dataDict['lat'][-1]]
        spatialplot(dataDict, bb=bb, timeidx = timeidx, title='vorticity')
    return dataDict

def divergence(dataDict,timeidx=0,plot = False):
    div = np.empty((len(dataDict['time']),len(dataDict['lat']), len(dataDict['lon'])))
    div[:] = np.nan
    for t in range(1, len(dataDict['time'])):
        for i in range(1,len(dataDict['lon'])-1):
            for j in range(1,len(dataDict['lat'])-1):
                dvdx = (dataDict['u'][t, j + 1, i + 1] - dataDict['u'][t, j - 1, i - 1]) / 2 * (dataDict['lon'][1] - dataDict['lon'][0])
                dudy = (dataDict['v'][t, j + 1, i + 1] - dataDict['v'][t, j - 1, i - 1]) / 2 * (dataDict['lat'][1] - dataDict['lat'][0])
                div[t,j, i] = dvdx + dudy
    dataDict['div'] = div
    if plot == True:
        bb = [dataDict['lon'][0], dataDict['lon'][-1], dataDict['lat'][0], dataDict['lat'][-1]]
        spatialplot(dataDict, bb=bb, timeidx = timeidx,title='divergence')
    return dataDict

def sshgrad(dataDict,timeidx=0,plot=False):
    sshgrad = np.empty((len(dataDict['timeReg']), len(dataDict['lat']), len(dataDict['lon'])))
    sshgrad[:] = np.nan
    for t in range(1, len(dataDict['timeReg'])):
        for i in range(1,len(dataDict['lon'])-1):
            for j in range(1,len(dataDict['lat'])-1):
                dhdx = (dataDict['ssh_ut'][t, j + 1, i + 1] - dataDict['ssh_ut'][t, j - 1, i - 1]) / 2 * (dataDict['lon'][1] - dataDict['lon'][0])
                dhdy = (dataDict['ssh_ut'][t, j + 1, i + 1] - dataDict['ssh_ut'][t, j - 1, i - 1]) / 2 * (dataDict['lat'][1] - dataDict['lat'][0])
                sshgrad[t,j, i] = dhdx + dhdy
    dataDict['sshgrad'] = sshgrad
    if plot == True:
        bb = [dataDict['lon'][0], dataDict['lon'][-1], dataDict['lat'][0], dataDict['lat'][-1]]
        spatialplot(dataDict, bb=bb, timeidx = timeidx,title='divergence')
    return dataDict

def residuals_utide(dataDict,timeidx=0):

    vort = np.empty((len(dataDict['timeReg']),len(dataDict['lat']), len(dataDict['lon'])))
    vort[:] = np.nan

    div = np.empty((len(dataDict['timeReg']),len(dataDict['lat']), len(dataDict['lon'])))
    div[:] = np.nan

    for t in range (1,len(dataDict['timeReg'])):
        for i in range(1,len(dataDict['lon'])-1):
            for j in range(1,len(dataDict['lat'])-1):
                dvdx = (dataDict['ut'][t, j + 1, i + 1] - dataDict['ut'][t, j - 1, i - 1]) / 2 * (dataDict['lon'][1] - dataDict['lon'][0])
                dudy = (dataDict['vt'][t, j + 1, i + 1] - dataDict['vt'][t, j - 1, i - 1]) / 2 * (dataDict['lat'][1] - dataDict['lat'][0])
                vort[t,j, i] = dvdx - dudy
                div[t, j, i] = dvdx + dudy
    dataDict['vort'] = vort
    dataDict['div'] = div
    return dataDict
