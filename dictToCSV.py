import pickle
import csv
import datetime as datetime
import numpy as np
from Tr_find import Tr_find
matdata = 'Tr_FORCE_2012-2031.mat'
fm = 'Aug142019.p'
infile = open(fm, 'rb')
data = pickle.load(infile,encoding='latin1')
infile.close()

stationStart = float(737651)
step = 738280.041666667-738280
stationEnd = float(737652)
timeStamp = np.linspace(stationStart, stationEnd, 24)


for i in range(24):
    data['time'][i,:,:] = np.tile(timeStamp[i], (30, 50))
data['Tr'] = np.zeros((24,30,50))
data['Tt'] = np.zeros((24,30,50))
stationStart = float(737651)
step = 738280.041666667-738280
stationEnd = float(737652)
timeVec = np.linspace(stationStart, stationEnd, 24)

for t in range(24):
    temp = Tr_find(matdata, timeVec[t])
    for i in range(30):
        for j in range(50):
            data['Tr'][t,i,j] = temp['Tr_val']
            data['Tt'][t, i, j] = temp['Tt_val']

for key in data.keys():
    data[key] = np.reshape(data[key],(36000))


with open('Aug142019.csv','w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(data.keys())
    writer.writerows(zip(*data.values()))

