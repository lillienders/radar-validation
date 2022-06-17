from get_data import *
import utide as ut
from val_plots import *
from calcs import *
from scipy import interpolate
from Tr_find import Tr_find
import scipy.io as sio
import seaborn as sns

# Path to files
radarFile = 'Data/data_at_dec18_adcp.p'
adcpFile = 'Data/adcp2018.mat'

# Load file w Radar and Model Data @ Dec 18 Deploy
infile = open(radarFile, 'rb')
dec_18_data = pickle.load(infile)
infile.close()

# Add times for harmonic regression
startTime = datetime.datetime(2017, 12, 14, 0, 0, 0)
step = datetime.datetime(2017, 12, 14, 0, 20, 0)
endTime = datetime.datetime(2018, 2, 22, 0, 0, 0)

stepTime = step - startTime

dec_18_data['time_adcp'] = date2num(np.arange(startTime, endTime, stepTime))

# Solve model and radar coefs
coef_model = ut.solve(dec_18_data['time'], dec_18_data['u_model'],dec_18_data['v_model'],
                                lat=45.363, conf_int='linear', epoch=datetime.datetime(1970,1,1))
coef_radar = ut.solve(dec_18_data['time'], dec_18_data['u_radar'],dec_18_data['v_radar'],
                                lat=45.363, conf_int='linear', epoch=datetime.datetime(1970,1,1))

ut_velocities_model = ut.reconstruct(dec_18_data['time_adcp'], coef_model, epoch=datetime.datetime(1970,1,1))
ut_velocities_radar = ut.reconstruct(dec_18_data['time_adcp'], coef_radar, epoch=datetime.datetime(1970,1,1))

dec_18_data['ut_model']=ut_velocities_model['u']
dec_18_data['vt_model']=ut_velocities_model['v']
dec_18_data['st_model']= np.sqrt(np.square(ut_velocities_model['u'])+np.square(ut_velocities_model['v']))

dec_18_data['ut_radar']=ut_velocities_radar['u']
dec_18_data['vt_radar']=ut_velocities_radar['v']
dec_18_data['st_radar']= np.sqrt(np.square(ut_velocities_radar['u'])+np.square(ut_velocities_radar['v']))

# Load adcp mat
adcpDict = getADCP(adcpFile,loc = [-64.42769444,45.363])
adcpDict['time_datetime'] = adcpDict['time'] - 719529

print(datetime.datetime.fromtimestamp(737043.83507))

fig, ax = plt.subplots(figsize=(20, 8))
cols = ['darkblue', 'crimson']
ax.plot(dec_18_data['time_adcp'], dec_18_data['st_radar'], color=cols[0], linewidth=3, zorder=0)
ax.plot(adcpDict['time_datetime'][25:6212], adcpDict['speed'][0:6187], color=cols[1], linewidth=3, zorder=1)

plt.xticks(fontsize=28, rotation=45)
plt.yticks(fontsize=28)
plt.gcf().subplots_adjust(bottom=0.25)
plt.ylabel('Speed (m/s)', fontsize=28)
ax.yaxis.set_label_coords(-0.08, 0.5)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%d-%b"))
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
ax.set_xlim(dec_18_data['time_adcp'][1300], dec_18_data['time_adcp'][1600])
ax.set_ylim(0,6.5)
ax.legend(['Radar', 'ADCP'],loc='upper center',fontsize=26)
plt.rcParams.update({'font.size': 28})
plt.savefig( 'todayisbrutal.png', format='png', bbox_inches="tight")
plt.show()

t_shift = adcpDict['time_datetime'][25:6212]
s_shift = adcpDict['speed'][0:6187]
u_shift = adcpDict['u'][0:6187]
v_shift = adcpDict['v'][0:6187]
f = interpolate.interp1d(t_shift, s_shift,fill_value="extrapolate")
f2 = interpolate.interp1d(t_shift, u_shift,fill_value="extrapolate")
f3 = interpolate.interp1d(t_shift, v_shift,fill_value="extrapolate")
dec_18_data['s_interp'] = f(dec_18_data['time_adcp'])
dec_18_data['u_interp'] = f2(dec_18_data['time_adcp'])
dec_18_data['v_interp'] = f3(dec_18_data['time_adcp'])

rho_num = (np.nansum(dec_18_data['u_interp'][100:4700].real * dec_18_data['ut_radar'][100:4700].real +
                          dec_18_data['v_interp'][100:4700].real * dec_18_data['vt_radar'][100:4700].real))
rho_denom = (np.sqrt(np.nansum((np.square(dec_18_data['u_interp'][100:4700].real) + np.square(dec_18_data['v_interp'][100:4700].real))))) * \
                    (np.sqrt(np.nansum((np.square(dec_18_data['ut_radar'][100:4700].real) + np.square(dec_18_data['vt_radar'][100:4700].real)))))
rho = rho_num/rho_denom
print(rho)
sum = np.sum(np.square(dec_18_data['st_radar'][100:4700]-dec_18_data['s_interp'][100:4700]))
RMSE_speed = np.sqrt((1/len(dec_18_data['time_adcp'][100:4700]))*sum)
print(RMSE_speed)