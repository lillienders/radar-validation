import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
from matplotlib.ticker import FormatStrFormatter
import numpy as np
"""
Functions to make radar validation plots

Provides functions:
    - plotTS: grab data from dictionary, return u-tide coefficients for each location
    - spatialPlot: grab u-tide coefficients and use them to reconstruct time-series at each location

Lilli Enders (lilli.enders@outlook.com)
June 2021
"""

def plotTS(data,vars,loc,xlim=None,title='',**kwargs):
    # Find index of input location in data
    lonidx = np.abs(data['lon'] - loc[0]).argmin()
    latidx = np.abs(data['lat'] - loc[1]).argmin()

    # Set up figure
    fig, ax = plt.subplots(figsize=(20,8))
    cols = ['darkblue','crimson','black','red']
    for i in range(len(vars)):
        if len(data[vars[i]]) == len(data['time']):
            ax.plot(data['time'], data[vars[i]][:,latidx,lonidx], color = cols[i], linewidth=3, zorder=i)
        else:
            ax.plot(data['timeReg'], data[vars[i]][:, latidx, lonidx], color = cols[i],linewidth=3, zorder=i)

    plt.xticks(fontsize=15, rotation=45)
    plt.gcf().subplots_adjust(bottom=0.25)
    plt.ylabel('Speed Residuals (m/s)', fontsize=20)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%d-%b"))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])
    plt.rcParams.update({'font.size': 22})
    plt.savefig( title+'.png', format='png', bbox_inches="tight")
    plt.show()
    return

def compareTS(radardata,modeldata,radarvar,modelvar,loc,xlim=None,title='',**kwargs):
    # Find index of input location in data
    lonidx = np.abs(radardata['lon'] - loc[0]).argmin()
    latidx = np.abs(radardata['lat'] - loc[1]).argmin()

    # Set up figure
    fig, ax = plt.subplots(figsize=(20,8))
    cols = ['darkblue','crimson','black','red']

    if len(radardata[radarvar]) == len(radardata['time']):
        ax.plot(radardata['time'], radardata[radarvar][:, latidx, lonidx], color=cols[0], linewidth=3, zorder=0)
    else:
        ax.plot(radardata['timeReg'], radardata[radarvar][:, latidx, lonidx], color=cols[0], linewidth=3, zorder=0)
    if len(modeldata[modelvar]) == len(modeldata['time']):
        ax.plot(modeldata['time'], modeldata[modelvar][:, latidx, lonidx], color=cols[1], linewidth=3, zorder=1)
    else:
        ax.plot(modeldata['timeReg'], modeldata[modelvar][:, latidx, lonidx], color=cols[1], linewidth=3, zorder=1)

    plt.xticks(fontsize=15, rotation=45)
    plt.gcf().subplots_adjust(bottom=0.25)
    plt.ylabel('Speed Residuals (m/s)', fontsize=20)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%d-%b"))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])
    plt.rcParams.update({'font.size': 22})
    plt.savefig( title+'.png', format='png', bbox_inches="tight")
    plt.show()
    return

def spatialplot(radardata,modeldata,title='',**kwargs):
    fig,ax = plt.subplots(2,figsize=(14,11))
    ax[0].set_facecolor('gainsboro')
    ax[1].set_facecolor('gainsboro')
    colorplot1 = ax[0].pcolor(modeldata['lon'],modeldata['lat'],modeldata['M2Grid'],vmax=5,vmin=0,shading='auto',cmap=plt.cm.jet)
    colorplot2 = ax[1].pcolor(radardata['lon'],radardata['lat'],radardata['M2Grid'],vmax=5,vmin=0,shading='auto',cmap=plt.cm.jet)
    ax[0].title.set_text('FVCOM')
    ax[1].title.set_text('Radar')
    plt.xlabel('Longitude (°)')
    ax[0].set_ylabel('Latitude (°)')
    ax[1].set_ylabel('Latitude (°)')
    fig.subplots_adjust(right=0.8)
    cbarax =fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(colorplot2,cax=cbarax)
    cbar.set_label('M2 Amplitude',rotation=90)
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 22})
    plt.savefig(title +'.png', dpi=1000)
    plt.show()
    return