import utide as ut
from get_data import *
from do_utide import *
from val_plots import *

import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
"""
Read in Files (Model, Radar) 
"""
modelFile = 'Data/FVCOMRadarGrid.p'
radarFile = 'Data/Radar150.p'

modelDict = getFVCOM(modelFile)
radarDict = getRadar(radarFile)

"""
Fit constituents (@ One Location) 
"""

radarCoef = ut.solve(radarDict['time'],radarDict['u'][:,1,1],radarDict['v'][:,1,1],
                lat=radarDict['lat'][1],conf_int='linear',epoch='python')
radarVel = ut.reconstruct(radarDict['timeReg'],radarCoef,epoch='python')

modelCoef = ut.solve(modelDict['time'],modelDict['u'][:,1,1],modelDict['v'][:,1,1],
                lat=radarDict['lat'][1],conf_int='linear',epoch='python')
modelVel = ut.reconstruct(modelDict['time'],modelCoef,epoch='python')

"""
Fit constituents (@ Multiple/All Locations on Radar Grid) 
"""
locsToGrab = [(-64.43059557,45.36149244),(-64.42868038,45.36284225),(-64.42868038,45.36284225)]
modelDict,modelCoefs = solveCoefs(modelDict)
radarDict,radarCoefs = solveCoefs(radarDict)

"""
Time Series Plots 
"""
'''
# Compare w/in one data source (ex compare raw data to reconstructed data)
plotTS(modelDict,vars = ['st','speed'],xlim=[modelDict['time'][500],modelDict['time'][800]],
       loc=[-64.43059557,45.36149244],title='test')

# Compare between data sources (ex radar data to model data)
compareTS(radarDict,modelDict,radarvar='st',modelvar = 'st',xlim=[modelDict['time'][500],modelDict['time'][800]],
          loc=[-64.43059557,45.36149244],title='')
'''
"""
Spatial Plots
"""
spatialplot(radarDict,modelDict,title='testM2')