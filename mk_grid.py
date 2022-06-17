from get_data import *
from do_utide import *
from calcs import *
import pandas as pd
from pylab import *

"""
Read in Files (Model, Radar) 
"""
modelFile = 'Data/fvcom_on_Radar_Grid_surface.p'
radarFile = 'Data/RadarMerged.p'

modelDict = getFVCOM(modelFile)
radarDict = getRadar(radarFile)

"""
Calculate Model Data
"""
modelDict,coefsModel = solveCoefs(modelDict,SSH=True)
print('Model Done')
modelDict = sshgrad(modelDict)

"""
Calculate Radar Data
"""
radarDict,coefsRadar = solveCoefs(radarDict)
print('Radar Done')
radarDict = residuals_utide(radarDict)
grid  = {}
grid['lat'] = modelDict['lat']
grid['lon'] = modelDict['lon']
grid['time'] = modelDict['timeReg']
grid['ssh'] = modelDict['ssh_ut']
grid['ssh_grad'] = modelDict['sshgrad']
grid['u_vel'] = radarDict['ut']
grid['v_vel'] = radarDict['vt']
grid['vort'] = radarDict['vort']
grid['div'] = radarDict['div']

print('Prepping File ....')
with open( 'Data/grid.p', 'wb') as handle:
    pickle.dump(grid, handle, protocol=pickle.HIGHEST_PROTOCOL)
print('Grid Dictionary Exported (.p)')

f = 'Data/grid.p'
infile = open(f, 'rb')
grid = pickle.load(infile,encoding='latin1')
infile.close()

lat = []
lon =[]
time = []
ssh = []
ssh_grad = []
u_vel = []
v_vel = []
vort = []
div = []

for i in range(1, len(grid['lon']) - 1):
    for j in range(1, len(grid['lat']) - 1):
        time_temp = grid['time']
        lat_temp = np.repeat(grid['lat'][j], len(grid['time']))
        lon_temp = np.repeat(grid['lon'][i], len(grid['time']))
        ssh_temp = grid['ssh'][:,j,i]
        ssh_grad_temp = grid['ssh_grad'][:, j, i]
        u_temp = grid['u_vel'][:, j, i]
        v_temp = grid['v_vel'][:, j, i]
        vort_temp = grid['vort'][:, j, i]
        div_temp = grid['div'][:, j, i]

        time = np.append(time, time_temp)
        lat = np.append(lat, lat_temp)
        lon = np.append(lon, lon_temp)
        ssh = np.append(ssh,ssh_temp)
        ssh_grad = np.append(ssh_grad, ssh_grad_temp)
        u_vel = np.append(u_vel, u_temp)
        v_vel = np.append(v_vel, v_temp)
        vort = np.append(vort, vort_temp)
        div = np.append(div, div_temp)

grid_concat = {}
grid_concat['lat'] = lat
grid_concat['lon'] = lon

grid_concat['time'] = time
grid_concat['ssh'] = ssh
grid_concat['ssh_grad'] = ssh_grad
grid_concat['u_vel'] = u_vel
grid_concat['v_vel'] = v_vel
grid_concat['vort'] = vort
grid_concat['div'] = div

print('Prepping File ....')
with open( 'Data/grid_concat.p', 'wb') as handle:
    pickle.dump(grid_concat, handle, protocol=pickle.HIGHEST_PROTOCOL)
print('Concatenated Samples Exported (.p)')

f = 'Data/grid_concat.p'
infile = open(f, 'rb')
grid_concat = pickle.load(infile,encoding='latin1')
infile.close()
grid_concat['time'] = num2date(grid_concat['time'])
grid_concat['speed'] = np.sqrt(np.square(grid_concat['u_vel'])+np.square(grid_concat['v_vel']))
print('Prepping File ....')
with open( 'Data/grid_concat2.p', 'wb') as handle:
    pickle.dump(grid_concat, handle, protocol=pickle.HIGHEST_PROTOCOL)
print('Concatenated Samples Exported (.p)')

for j in range(len(grid_concat['time'])):
    grid_concat['time'][j] = grid_concat['time'][j].strftime("%m/%d/%Y, %H:%M")

griddf = pd.DataFrame.from_dict(grid_concat)
compression_opts = dict(method='zip',
                        archive_name='grid.csv')
griddf.to_csv('grid.zip', index=False,
          compression=compression_opts)
print('CSV Export Complete')
print('Done')
