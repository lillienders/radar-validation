from get_data import *
from calcs import *
from do_utide import *
import statistics as stats
import numpy as np
import math as math

"""
Error metrics used for radar validation. Calculates RMSE & NRMSE of speed metric, complex correlation, and phase angles

Lilli Enders (lilli.enders@outlook.com)
June 2022
"""
modelFile = 'Data/fvcom_on_Radar_Grid_surface.p'
radarFile = 'Data/RadarMerged.p'

modelDict = getFVCOM(modelFile)
radarDict = getRadar(radarFile)

modelDict,modelCoefs = solveCoefs(modelDict)
radarDict,radarCoefs = solveCoefs(radarDict)

RMSE_speed = np.zeros((len(radarDict['lat']),len(radarDict['lon'])))
n = np.zeros((len(radarDict['lat']),len(radarDict['lon'])))
mean_speed = np.zeros((len(radarDict['lat']),len(radarDict['lon'])))
rho = np.zeros((len(radarDict['lat']),len(radarDict['lon'])))
alpha = np.zeros((len(radarDict['lat']),len(radarDict['lon'])))

for i in range(len(radarDict['lat'])):
    for j in range(len(radarDict['lon'])):
        n[i,j] = sum(~np.isnan(radarDict['speed'][:, i, j]))
        mean_speed[i,j] = stats.mean(radarDict['speed'][:,i,j])

        rho_num = (np.nansum(modelDict['u'][:, i, j].real * radarDict['u'][:, i, j].real +
                          modelDict['v'][:, i, j].real * radarDict['v'][:, i, j].real))
        rho_denom = (np.sqrt(
            np.nansum((np.square(modelDict['u'][:, i, j].real) + np.square(modelDict['v'][:, i, j].real))))) * \
                    (np.sqrt(
            np.nansum((np.square(radarDict['u'][:, i, j].real) + np.square(radarDict['v'][:, i, j].real)))))
        rho[i,j] = rho_num / rho_denom

        alpha[i,j] = np.nansum(modelDict['u'][:, i, j].real * radarDict['v'][:, i, j].real -
                               modelDict['v'][:, i, j].real * radarDict['u'][:, i, j].real)/\
                     np.nansum(modelDict['u'][:, i, j].real * radarDict['u'][:, i, j].real +
                               modelDict['v'][:, i, j].real*radarDict['v'][:, i, j].real)

sum = np.zeros((len(radarDict['lat']),len(radarDict['lon'])))

for i in range(30):
    for j in range(50):
        alpha[i,j] = math.degrees(alpha[i,j])

radarDict['speed'] = np.nan_to_num(radarDict['speed'])
radarDict['u'] = np.nan_to_num(radarDict['u'])
radarDict['v'] = np.nan_to_num(radarDict['v'])
modelDict['speed'] = np.nan_to_num(modelDict['speed'])
modelDict['u'] = np.nan_to_num(modelDict['u'])
modelDict['v'] = np.nan_to_num(modelDict['v'])

for t in range(len(radarDict['time'])):
    sum = sum + np.square(radarDict['speed'][t,:,:]-modelDict['speed'][t,:,:])

RMSE_speed = np.sqrt((1/n)*sum)

for i in range(len(radarDict['lat'])):
    for j in range(len(radarDict['lon'])):
        if RMSE_speed[i,j] == np.Inf:
            RMSE_speed[i, j] = np.nan

NRMSE_speed = RMSE_speed/mean_speed
radarDict['RMSE_speed'] = RMSE_speed
radarDict['NRMSE_speed'] = NRMSE_speed
radarDict['complex_corr'] = rho
radarDict['phase_angle'] = alpha

bb = [-64.50,-64.35,45.32,45.38]
spatialplot(radarDict,bb,title='complex_corr',timeidx = 0)