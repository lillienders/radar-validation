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
bb = [-64.445,np.max(radarDict['lon']),45.35,np.max(radarDict['lat'])]
plot = spatialcompare(radarDict,modelDict,var = 'speed',bb=bb,timeidx=120)

date_cols = ["deploy_date_time","recover_date_time"]
df = pd.read_csv ('Data/MP_Stations.csv',parse_dates=date_cols)

"""
Set up Samples, loop through each row and assign variables to dictionary of dictionaries
"""
samples = {}
for index, row in df.iterrows():
    startTime = row['deploy_date_time']
    endTime = row['recover_date_time']
    step = datetime.timedelta(hours=1)
    radarDict['timeReg'] = date2num(np.arange(startTime, endTime, step))
    modelDict['timeReg'] = date2num(np.arange(startTime, endTime, step))
    locs = [(row['mean_long']-0.00195,row['mean_lat']),(row['mean_long']+0.00195,row['mean_lat']),
            (row['mean_long'],row['mean_lat']-0.00135),(row['mean_long'],row['mean_lat']+0.00135),
            (row['mean_long'],row['mean_lat'])]

    lonidx = np.abs(modelDict['lon'] - row['mean_long']).argmin()
    latidx = np.abs(modelDict['lat'] - row['mean_lat']).argmin()
    modelDict, coefsModel = solveCoefs(modelDict, SSH=True,locs = locs)
    modelDict = sshgrad(modelDict)

    radarDict, coefsRadar = solveCoefs(radarDict,locs = locs)

    samples[row['station']] = {}
    samples[row['station']]['lat'] = row['mean_lat']
    samples[row['station']]['lon'] = row['mean_long']
    samples[row['station']]['time'] = date2num(np.arange(startTime, endTime, step))
    samples[row['station']]['ssh'] = modelDict['ssh_ut'][:,latidx, lonidx]
    samples[row['station']]['ssh_grad'] = modelDict['sshgrad'][:,latidx-1, lonidx-1]
    samples[row['station']]['u_vel'] = radarDict['ut'][:,latidx, lonidx]
    samples[row['station']]['v_vel'] = radarDict['vt'][:,latidx, lonidx]
    samples[row['station']]['vort'], samples[row['station']]['div'] = residuals_utide(radarDict, lon=row['mean_long'], lat=row['mean_lat'])
    print(row['station'])

"""
Optional: Export samples dictionary, so that we don't have to repeat computation later
"""
print('Prepping File ....')
with open( 'Data/samples.p', 'wb') as handle:
    pickle.dump(samples, handle, protocol=pickle.HIGHEST_PROTOCOL)
print('Sample Dictionary Exported (.p)')

"""
Concatenate vars to export to csv
"""
station = []
lat = []
lon =[]

time = []
ssh = []
ssh_grad = []
u_vel = []
v_vel = []
vort = []
div = []

for index, row in df.iterrows():
    station_temp = np.repeat([row['station']],len(samples[row['station']]['time']))
    lat_temp = np.repeat(samples[row['station']]['lat'],len(samples[row['station']]['time']))
    lon_temp = np.repeat(samples[row['station']]['lon'], len(samples[row['station']]['time']))

    station = np.append(station,station_temp)
    lat = np.append(lat,lat_temp)
    lon = np.append(lon,lon_temp)

    time = np.append(time,samples[row['station']]['time'])
    ssh = np.append(ssh, samples[row['station']]['ssh'])
    ssh_grad = np.append(ssh_grad, samples[row['station']]['ssh_grad'])
    u_vel = np.append(u_vel, samples[row['station']]['u_vel'])
    v_vel = np.append(v_vel, samples[row['station']]['v_vel'])
    vort = np.append(vort, samples[row['station']]['vort'])
    div = np.append(div, samples[row['station']]['div'])

samples_concat = {}
samples_concat['station'] = station
samples_concat['lat'] = lat
samples_concat['lon'] = lon

samples_concat['time'] = time
samples_concat['ssh'] = ssh
samples_concat['ssh_grad'] = ssh_grad
samples_concat['u_vel'] = u_vel
samples_concat['v_vel'] = v_vel
samples_concat['vort'] = vort
samples_concat['div'] = div

"""
Dump concatenated files (just in case you need it)
"""
print('Prepping File ....')
with open( 'Data/samples_concat.p', 'wb') as handle:
    pickle.dump(samples_concat, handle, protocol=pickle.HIGHEST_PROTOCOL)
print('Concatenated Samples Exported (.p)')

# Convert time to string (optional - otherwise python datetime format)
samples_concat['time'] = num2date(samples_concat['time'])

for j in range(len(samples_concat['time'])):
    samples_concat['time'][j] = samples_concat['time'][j].strftime("%m/%d/%Y, %H:%M")

samplesdf = pd.DataFrame.from_dict(samples_concat)
compression_opts = dict(method='zip',
                        archive_name='samples.csv')
samplesdf.to_csv('samples.zip', index=False,
          compression=compression_opts)
print('CSV Export Complete')

