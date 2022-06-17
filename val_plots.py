from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mtick
from matplotlib.ticker import FormatStrFormatter
import matplotlib.tri as Tri
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
import seaborn as sns
import matplotlib.gridspec as gridspec
import numpy as np
import pickle
import string

"""
Functions to make radar validation plots
Provides functions:
    - hex_to_rgb: Converts hex codes to rgb for plotting w pyplot 
    - rgb_to_dec: Converts rgb to decimal for plotting w pyplot
    - get_continuous_cmap: creates custom colourmap given input colours in dec
    - plotTS: Plots time series of input variable(s) at one location for radar OR model 
    - compareTS: Plots two time series: one of radar and one of model given input variable(s) 
    - spatialplot: Creates spatial plot of input parameter for radar or model data over whole grid 
    - radarcount: Creates spatial plot of number of radar measurements that pass ellipse QC test 
    - plotellipse: Plots values of ellipse parameter (error measure) given input radar data
Lilli Enders (lilli.enders@outlook.com)
June 2022
"""
def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def plotTS(data,vars,loc,xlim=None,title='',**kwargs):
    # Find index of input location in data
    lonidx = np.abs(data['lon'] - loc[0]).argmin()
    latidx = np.abs(data['lat'] - loc[1]).argmin()

    # Set up figure
    fig, ax = plt.subplots(figsize=(15,6))
    cols = ['midnightblue','crimson','red']
    for i in range(len(vars)):
        if len(data[vars[i]]) == len(data['time']):
            ax.scatter(data['time'], data[vars[i]][:,latidx,lonidx], color = cols[i], zorder=i)
        else:
            ax.plot(data['timeReg'], data[vars[i]][:, latidx, lonidx], color = cols[i],linewidth=1, zorder=i)

    plt.xticks(fontsize=20, rotation=45)
    plt.yticks(fontsize=20)
    plt.gcf().subplots_adjust(bottom=0.25)
    plt.ylabel('Signed Speed (m/s)', fontsize=20)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%d-%b"))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    #ax.legend(['Reconstructed TS', 'Raw TS'],loc='upper center',fontsize=16)
    if xlim != None:
        ax.set_xlim(xlim[0], xlim[1])
    ax.set_xlim(data['timeReg'][0], data['timeReg'][-1])
    ax.set_ylim(0, 5.5)
    plt.rcParams.update({'font.size': 20})
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


def spatialplot(data,bb,title='',timeidx = 0,**kwargs):
    infile = open('Data/tryp.pickle', 'rb')
    bathydict = pickle.load(infile)
    infile.close()
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_facecolor('white')
    #ax.set_xlim(bb[0] - 0.005, bb[1] + 0.005)
    #ax.set_ylim(bb[2] - 0.005, bb[3] + 0.005)
    ax.set_xlim(bb[0], bb[1])
    ax.set_ylim(bb[2], bb[3])
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.xlabel('Longitude (°)')
    plt.ylabel('Latitude (°)')
    ax.xaxis.set_label_coords(0.5, -0.075)
    ax.yaxis.set_label_coords(-0.075, 0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tripcolor(tri, bathydict['h'], cmap=plt.cm.Greys, vmin=-500, vmax=2000)
    #colorplot = ax.pcolormesh(data['lon'], data['lat'], data['vort'][timeidx,:,:],
    #                          vmin=-0.0025, vmax=0.0025,
                              #cmap=plt.cm.get_cmap("gist_earth"), shading='auto', zorder=10)
    #                          cmap = plt.cm.get_cmap("bwr"), shading = 'auto', zorder = 10)

    hex_list = ['#2c5a8c', '#FFFFFF', '#cf3f47']
    #colorplot = ax.pcolormesh(data['lon'], data['lat'], data['diff'][timeidx,:,:],
    colorplot=ax.pcolormesh(data['lon'], data['lat'], data['M2thetadiff'],
                              vmin=-5, vmax=5,
                              #cmap=plt.cm.get_cmap("gist_earth"), shading='auto', zorder=10)
                              #cmap = plt.cm.get_cmap("Spectral_r"), shading ='auto', zorder = 10)
                              cmap=get_continuous_cmap(hex_list),shading ='auto', zorder = 10)

    #colorplot = ax.contourf(data['lon'], data['lat'], data['complex_corr'], 100,
    #                          vmin=0, vmax=1,
    #                          cmap=plt.cm.get_cmap("Spectral_r"), zorder=10)
                             #cmap = get_continuous_cmap(hex_list), zorder = 10)
    for item in ([ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(20)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    fig.subplots_adjust(right=0.82, bottom=0.15)
    cbarax = fig.add_axes([0.875, 0.125, 0.025, 0.75])
    cbar = plt.colorbar(colorplot, cax=cbarax)
    cbar.outline.set_visible(False)
    cbar.set_label('bathy std (m)', rotation=270, size='20', labelpad=25)
    cbar.ax.tick_params(labelsize='20')
    CLA_loc = [-64.44090, 45.36150]  # At south-west corner
    CLA = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899,
                        fill=False, color='k', linewidth=3.0,zorder=30)
    ax.add_patch(CLA)
    '''
    ADCP_dict = {'2013_600_1': [-64.42426667, 45.36698333], '2013_400_1': [-64.42543333, 45.36653333],
                 '2012_600_8': [-64.43078333, 45.36486667], '2012_400_2': [-64.42926667, 45.36508333],
                 '2012_600_7': [-64.41861667, 45.3658], '2012_600_6': [-64.41196667, 45.36826667],
                 '2012_600_5': [-64.42025, 45.36235], '2012_600_4': [-64.408, 45.36506667],
                 '2012_600_3': [-64.43736667, 45.36558333], '2012_400_1': [-64.43563333, 45.36535],
                 '2012_600_2': [-64.41513333, 45.3671], '2012_600_1': [-64.40541667, 45.36755],
                 '2011_400_2': [-64.42435, 45.36513333], '2011_600_2': [-64.42193333, 45.36428333],
                 '2011_400_1': [-64.43735, 45.36513333], '2011_600_1': [-64.4299, 45.36478333],
                 '2008_4': [-64.44388889, 45.36972222], '2008_3': [-64.44805556, 45.37083333],
                 '2008_2': [-64.46388889, 45.35611111], '2008_1': [-64.45888889, 45.36472222],
                 '2015_400_1': [-64.4047222, 45.369444], '2015_400_2': [-64.405, 45.369444],
                 '2015_400_3': [-64.405, 45.3697222], '2018_500_1': [-64.42769444, 45.363],
                 '2018_500_2': [-64.42747222, 45.36313889], '2018_500_3': [-64.42165, 45.364467],
                 '2018_500_4': [-64.422, 45.36375], '2018_500_6': [-64.42775, 45.36319444]}
    plot_adcp = {"lon": [], "lat": [], "label": []}
    for label, i in ADCP_dict.items():
        plot_adcp["lon"].append(i[0])
        plot_adcp["lat"].append(i[1])
        plot_adcp["label"].append(label)
    ax.scatter(plot_adcp["lon"], plot_adcp["lat"], color='k', marker='.', s=100, zorder=15)
    '''
    #levels = [20,40,60]
    #isolines = ax.tricontour(tri, bathydict['h'], colors='w', linewidths=2, levels=levels, zorder=20)
    #plt.clabel(isolines, levels=[20,40,60], fontsize='14', fmt='%.0f', inline=True,zorder=30)
    #for i in range(1, len(data['lon']) - 1):
    #    for j in range(1, len(data['lat']) - 1):
    #        ax.quiver(data['lon'][i], data['lat'][j], data['v'][timeidx, j, i],
    #                  data['v'][timeidx, j, i], color='k', headwidth=2, width=0.0015, scale=60, zorder=20)
    #fig.suptitle(mdates.DateFormatter("%Y-%m-%d %H:%M")(data['time'][timeidx]),
    #            fontsize='16',y=0.92)
    #plt.savefig('vorticity_animation/vorticity'+str(timeidx)+'.png', dpi=300)
    #plt.savefig('vorticitybwr_short/vorticity' + str(timeidx) + '.png', dpi=300)
    #plt.close()
    plt.savefig(title+'.png', dpi=1000)
    plt.show()

def spatialcompare(radardata,modeldata,var,bb,timeidx=0,title='',**kwargs):

    infile = open('Data/tryp.pickle', 'rb')
    bathydict = pickle.load(infile)
    infile.close()
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)

    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    fig, ax = plt.subplots(1,2,figsize=(32,10))

    ax[0].set_facecolor('white')
    ax[0].set_xlim(bb[0], bb[1])
    ax[0].set_ylim(bb[2], bb[3])
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    ax[1].set_facecolor('white')
    #ax[1].set_xlim(bb[0] - 0.005, bb[1] + 0.005)
    #ax[1].set_ylim(bb[2] - 0.005, bb[3] + 0.005)
    ax[1].set_xlim(bb[0], bb[1])
    ax[1].set_ylim(bb[2], bb[3])
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].set(xlabel='Longitude (°)',ylabel='Latitude (°)')
    ax[1].set(xlabel='Longitude (°)', ylabel='Latitude (°)')
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[0].tripcolor(tri, bathydict['h'], cmap=plt.cm.Greys, vmin=-500, vmax=2000)
    ax[1].tripcolor(tri, bathydict['h'], cmap=plt.cm.Greys, vmin=-500, vmax=2000)
    colorplot1 = ax[0].contourf(modeldata['lon'], modeldata['lat'], modeldata[var], 100,
                              vmin=0, vmax=5,
                              cmap=plt.cm.get_cmap("Spectral_r"), zorder=10)
    colorplot2 = ax[1].contourf(modeldata['lon'], modeldata['lat'], radardata[var], 100,
                              vmin=0, vmax=5, cmap=plt.cm.get_cmap("Spectral_r"), zorder=10)
    #colorplot1 = ax[0].pcolor(modeldata['lon'], modeldata['lat'], modeldata[var][timeidx+243,:,:],
    #                          vmin=-5, vmax=5,shading='auto', cmap=plt.cm.get_cmap("bwr"))

    #colorplot2 = ax[1].pcolor(radardata['lon'], radardata['lat'], radardata[var][timeidx,:,:],
    #                           vmin=-5, vmax=5,shading='auto', cmap=plt.cm.get_cmap("bwr"))
    for item in ([ax[0].xaxis.label, ax[0].yaxis.label,ax[1].xaxis.label, ax[1].yaxis.label]):
        item.set_fontsize(22)
    for item in (ax[0].get_xticklabels() + ax[0].get_yticklabels() + ax[1].get_xticklabels() + ax[1].get_yticklabels()):
        item.set_fontsize(20)
    '''
    levels = [15,40]
    isolines1 = ax[0].tricontour(tri, bathydict['h'], colors='k', linewidths=3, levels=levels, zorder=20)
    isolines2 = ax[1].tricontour(tri, bathydict['h'], colors='k', linewidths=3, levels=levels, zorder=20)
    ax[0].clabel(isolines1, levels=[15,40], fontsize='20', fmt='%.0f', inline=True,zorder=30)
    ax[1].clabel(isolines2, levels=[15,40], fontsize='20', fmt='%.0f', inline=True, zorder=30)
    '''
    CLA_loc = [-64.44090, 45.36150]  # At south-west corner
    CLA = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899,
                        fill=False, color='k', linewidth=2.0, zorder=10)
    CLA2 = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899,
                        fill=False, color='k', linewidth=2.0, zorder=10)
    ax[0].add_patch(CLA)
    ax[1].add_patch(CLA2)

    '''
    for i in range(1, len(modeldata['lon']) - 1):
        for j in range(1, len(modeldata['lat']) - 1):
            
            ax[0].quiver(modeldata['lon'][i], modeldata['lat'][j], modeldata['u'][timeidx + 243, j, i],
                      modeldata['v'][timeidx, j, i], color='k', zorder=20)
            ax[1].quiver(radardata['lon'][i], radardata['lat'][j], radardata['u'][timeidx, j, i],
                      radardata['v'][timeidx, j, i], color='k', zorder=20)
           
            ax[0].quiver(modeldata['lon'][i], modeldata['lat'][j], modeldata['u'][timeidx + 243, j, i],
                      modeldata['v'][timeidx, j, i], color='k', headwidth=2, width=0.0015, scale=50, zorder=20)
            ax[1].quiver(radardata['lon'][i], radardata['lat'][j], radardata['u'][timeidx, j, i],
                      radardata['v'][timeidx, j, i], color='k', headwidth=2, width=0.0015, scale=50, zorder=20)
    '''

    for n, ax in enumerate(fig.get_axes()):
        ax.text(-0.1, 1.01, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,
                size=30,weight='bold',style='italic')

    fig.subplots_adjust(right=0.90, bottom=0.15)
    cbarax = fig.add_axes([0.925, 0.125, 0.025, 0.75])
    #cbarax = fig.add_axes()

    cbar = plt.colorbar(colorplot1, cax=cbarax)
    cbar.outline.set_visible(False)
    cbar.set_label('M2 Amplitude (m)', rotation=270, size='22', labelpad=25)
    cbar.ax.tick_params(labelsize='20')

    plt.subplots_adjust(left=0.05,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.125,
                        hspace=0.35)
    plt.savefig('M2Comp.png', dpi=300)
    plt.show()
    return


def spatialcompare2(data, vars, bb, timeidx=0, title='', **kwargs):
    infile = open('Data/tryp.pickle', 'rb')
    bathydict = pickle.load(infile)
    infile.close()
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)

    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    fig, ax = plt.subplots(1, 2, figsize=(20,6))

    ax[0].set_facecolor('white')
    ax[0].set_xlim(bb[0], bb[1])
    ax[0].set_ylim(bb[2], bb[3])
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].set_facecolor('white')
    ax[1].set_xlim(bb[0], bb[1])
    ax[1].set_ylim(bb[2], bb[3])
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].set(xlabel='Longitude (°)', ylabel='Latitude (°)')
    ax[1].set(xlabel='Longitude (°)', ylabel='Latitude (°)')
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[0].tripcolor(tri, bathydict['h'], cmap=plt.cm.Greys, vmin=-500, vmax=2000)
    ax[1].tripcolor(tri, bathydict['h'], cmap=plt.cm.Greys, vmin=-500, vmax=2000)

    hex_list = ['#2c5a8c', '#FFFFFF', '#cf3f47']
    colorplot1 = ax[0].pcolor(data['lon'], data['lat'], data[vars[0]]*100,
                                vmin=-50, vmax=50,
                                #cmap=plt.cm.get_cmap("Spectral_r"), zorder=10)
                                cmap=get_continuous_cmap(hex_list), shading = 'auto', zorder=10)


    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = fig.colorbar(colorplot1, cax=cax, orientation='vertical')
    cbar.outline.set_visible(False)
    #cbar.set_label('M2 Amplitude Diff (m)', rotation=270, size='22', labelpad=25)
    cbar.ax.tick_params(labelsize='20')
    cbar.ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))


    colorplot2 = ax[1].pcolor(data['lon'], data['lat'], data[vars[1]],
                               vmin=-10, vmax= 10,
                              cmap=get_continuous_cmap(hex_list), shading='auto', zorder=10)
                              #cmap=plt.cm.get_cmap("Spectral_r"), zorder=10)


    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = fig.colorbar(colorplot2, cax=cax, orientation='vertical')
    cbar.outline.set_visible(False)
    #cbar.set_label('M2 Phase Diff (°)', rotation=270, size='22', labelpad=25)
    cbar.ax.tick_params(labelsize='20')
    cbar.ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))

    for item in ([ax[0].xaxis.label, ax[0].yaxis.label, ax[1].xaxis.label, ax[1].yaxis.label]):
        item.set_fontsize(20)
    for item in (ax[0].get_xticklabels() + ax[0].get_yticklabels() + ax[1].get_xticklabels() + ax[1].get_yticklabels()):
        item.set_fontsize(20)

    CLA_loc = [-64.44090, 45.36150]  # At south-west corner
    CLA = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899,
                        fill=False, color='k', linewidth=2.0, zorder=10)
    CLA2 = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899,
                         fill=False, color='k', linewidth=2.0, zorder=10)
    ax[0].add_patch(CLA)
    ax[1].add_patch(CLA2)

    ''' 
    fig.subplots_adjust(right=0.90, bottom=0.15)
    cbarax = fig.add_axes([0.925, 0.125, 0.025, 0.75])
    # cbarax = fig.add_axes()

    cbar = plt.colorbar(colorplot1, cax=cbarax)
    cbar.outline.set_visible(False)
    cbar.set_label('M2 Amplitude Diff (m)', rotation=270, size='22', labelpad=25)
    cbar.ax.tick_params(labelsize='20')
        
    '''
    fig.tight_layout()

    plt.savefig('M2Comp.png', dpi=300)
    plt.show()
    return

def radarcount(data,title='',**kwargs):
    pointCount = np.empty((len(data['lat']), len(data['lon'])))
    pointCount[:] = np.nan
    for i in range(len(data['lon'])):
        for j in range(len(data['lat'])):
            uTemp = data['v'][:,j,i]
            #uTemp = data['M2maj'][j, i]
            if np.count_nonzero(~np.isnan(uTemp)) == 0:
                pointCount[j, i] = np.nan
            else:
                pointCount[j, i] = np.count_nonzero(~np.isnan(uTemp))

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_facecolor('gainsboro')
    colorplot = ax.pcolor(data['lon'], data['lat'], pointCount,shading='auto',
                          cmap=plt.cm.jet)
    #ax.title.set_text('FVCOM')
    plt.xlabel('Longitude (°)')
    ax.set_ylabel('Latitude (°)')
    fig.subplots_adjust(right=0.8)
    cbarax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(colorplot, cax=cbarax)
    cbar.set_label('No. Data Points', rotation=90)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.rcParams.update({'font.size': 22})
    plt.savefig(title + '.png', dpi=1000)

    plt.show()

def plotellipse(radardata,modeldata):
    plt.figure(figsize=(14, 8))
    ax = plt.gca()
    ax.set_facecolor('gainsboro')
    colorplot = ax.pcolor(modeldata['lon'], modeldata['lat'], modeldata['speed'][10, :, :], shading='auto',
                          cmap=plt.cm.jet)

    for i in range(len(radardata['lon'])):
        for j in range(len(radardata['lat'])):
            ellipseradar = Ellipse(xy=(radardata['lon'][i], radardata['lat'][j]), width=0.0005,
                                   height=radardata['M2maj'][j, i] / 1500, angle=radardata['M2phase'][j, i],
                                   edgecolor='k', fc='None', lw=2)
            ellipsemodel = Ellipse(xy=(modeldata['lon'][i], modeldata['lat'][j]), width=0.0005,
                                   height=modeldata['M2maj'][j, i] / 1500, angle=modeldata['M2phase'][j, i],
                                   edgecolor='r', fc='None', lw=2)
            ellipseADCP = Ellipse(xy=(-64.42747222, 45.363), width=0.0005, height=3.420 / 1500, angle=154.630,
                                  edgecolor='darkblue', fc='None', lw=2)
            ax.add_patch(ellipseradar)
            ax.add_patch(ellipsemodel)
            ax.add_patch(ellipseADCP)
    # ax.set_xlim(np.min(radarDict['lon']),np.max(radarDict['lon']))
    # ax.set_ylim(np.min(radarDict['lat']),np.max(radarDict['lat']))
    ax.set_xlim(-64.435, -64.425)
    ax.set_ylim(45.36, 45.365)
    plt.savefig('PreliminaryEllipse.png', dpi=1000)
    plt.show()


# OTHER PLOTS: Which I used for niche things that probably won't be repeated 
def M2quad(radardata,modeldata,bb,title='',**kwargs):
    infile = open('tryp.pickle', 'rb')
    bathydict = pickle.load(infile)
    infile.close()
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)

    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2,figsize = (14,8))
    cmap = plt.cm.get_cmap("gist_earth")
    #cmap = plt.cm.get_cmap("jet")
    cptl = ax1.pcolormesh(modeldata['lon'], modeldata['lat'], modeldata['M2maj'],
                               vmin=0, vmax=5,
                               cmap=cmap, shading='auto', zorder=10)
    cptr = ax2.pcolormesh(radardata['lon'], radardata['lat'], radardata['M2maj'],
                               vmin=0, vmax=5,
                               cmap=cmap, shading='auto', zorder=10)
    cpbl = ax3.pcolormesh(modeldata['lon'], modeldata['lat'], modeldata['M2phase'],
                               vmin=100, vmax=170,
                               cmap=cmap, shading='auto', zorder=10)
    cpbr = ax4.pcolormesh(radardata['lon'], radardata['lat'], radardata['M2phase'],
                               vmin=100, vmax=170,
                               cmap=cmap, shading='auto', zorder=10)
    for ax in fig.get_axes():
        ax.set_facecolor('white')
        ax.set_xlim(bb[0] + 0.001, bb[1] - 0.001)
        ax.set_ylim(bb[2], bb[3])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        #ax.xaxis.set_label_coords(0.5, -0.075)
        #ax.yaxis.set_label_coords(-0.075, 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tripcolor(tri, bathydict['h'], cmap=plt.cm.Greys, vmin=-500, vmax=2000)
        for item in ([ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(14)
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(12)
        isolines = ax.tricontour(tri, bathydict['h'], colors='w', linewidths=3,
                                       levels=[20, 40], zorder=20)
        plt.clabel(isolines, levels=[20, 40], colors = 'k',fontsize='12', fmt='%.0f', inline=True)

    ax1.get_xaxis().set_ticks([])
    ax2.get_xaxis().set_ticks([])

    for n, ax in enumerate(fig.get_axes()):
        ax.text(-0.1, 1.01, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,
                size=16,weight='bold',style='italic')

    fig.subplots_adjust(right=0.84, bottom=0.15)
    cax1 = fig.add_axes([0.875, 0.55, 0.025, 0.32])
    cax2 = fig.add_axes([0.875, 0.15, 0.025, 0.32])
    cb1 = plt.colorbar(cptl, cax=cax1)
    cb2 = plt.colorbar(cpbl, cax=cax2)
    cb1.outline.set_visible(False)
    cb2.outline.set_visible(False)



    #cbar.set_label('M2 Amplitude (m/s)', rotation=270, size='14', labelpad=25)
    #cbar.ax.tick_params(labelsize='12')
    CLA_loc = [-64.44090, 45.36150]  # At south-west corner
    CLA = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899, fill=False, color='w', linewidth=2.0, zorder=10)
    for ax in fig.get_axes():
        CLA = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899, fill=False, color='w', linewidth=2.0, zorder=10)
        ax.add_patch(CLA)
    plt.savefig(title+'.png', dpi=1000)
    plt.show()

def spatialplotTS(data,bb,title='',timeidx = 0,**kwargs):
    fig, (ax0, ax1,ax2) = plt.subplots(3, 1, figsize=(16,14), gridspec_kw={'height_ratios': [5, 1,1]})
    infile = open('Data/tryp.pickle', 'rb')
    bathydict = pickle.load(infile)
    infile.close()
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)
    nodes = np.transpose(bathydict['nv'].data) - 1
    tri = Tri.Triangulation(bathydict['lon'], bathydict['lat'], triangles=nodes)
    ax0.set_facecolor('white')
    ax0.set_xlim(bb[0], bb[1])
    ax0.set_ylim(bb[2], bb[3])
    ax0.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax0.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax0.set(xlabel='Longitude (°)', ylabel = 'Latitude (°)')
    ax0.xaxis.set_label_coords(0.5, -0.075)
    ax0.yaxis.set_label_coords(-0.075, 0.5)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.tripcolor(tri, bathydict['h'], cmap=plt.cm.Greys, vmin=-500, vmax=2000)
    colorplot = ax0.pcolormesh(data['lon'], data['lat'], data['M2diff'],
                              #vmin=-3, vmax=3,
                              cmap=plt.cm.get_cmap("Spectral_r"),shading='auto', zorder=10)
    for item in ([ax0.xaxis.label, ax0.yaxis.label]):
        item.set_fontsize(16)
    for item in (ax0.get_xticklabels() + ax0.get_yticklabels()):
        item.set_fontsize(14)
    fig.subplots_adjust(hspace = 0.4)
    #cbarax = ax0.add_axes([0.875, 0.125, 0.025, 0.75])
    cbarax = ax0.inset_axes([1.01, 0,0.03, 1], transform=ax0.transAxes)
    cbar = fig.colorbar(colorplot, ax=ax0, cax=cbarax)
    #cbar = ax0.colorbar(colorplot, cax=ax0)
    cbar.outline.set_visible(False)
    cbar.set_label('Along Stream Velocity', rotation=270, size='14', labelpad=25)
    cbar.ax.tick_params(labelsize='12')
    CLA_loc = [-64.44090, 45.36150]  # At south-west corner
    CLA = plt.Rectangle((CLA_loc[0], CLA_loc[1]), 0.02047, 0.00899,
                        fill=False, color='k', linewidth=3.0,zorder=30)
    ax0.add_patch(CLA)
    ''' 
    for i in range(1, len(data['lon']) - 1):
        for j in range(1, len(data['lat']) - 1):
            ax0.quiver(data['lon'][i], data['lat'][j], data['u'][timeidx, j, i],
                      data['v'][timeidx, j, i], color='k', headwidth=2, width=0.0015, scale=70, zorder=20)
    '''
    #c1 = hex_to_rgb(['#dd4c46'])
    ax1.plot(data['time'][250:400], data['u'][250:400, 15,15], color = '#dd4c46',linewidth=3)
    ax1.plot(data['time'][timeidx], data['u'][timeidx,15,15], "ro",color='black')
    ax1.set_xlim(data['time'][250], data['time'][399])
    ax1.set_title('Along Stream Velocity (m/s)',size=16,fontweight='bold')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%d-%b"))
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    ax1.set(xlabel='Time Stamp', ylabel='Velocity (m/s)')
    ax2.plot(data['time'][250:400], data['SSH'][250:400, 15,15],color = '#f88b52', linewidth=3)
    ax2.plot(data['time'][timeidx], data['SSH'][timeidx, 15, 15],"ro",color='black')
    ax2.set_xlim(data['time'][250], data['time'][399])
    ax2.set_title('Sea Surface Height (m)',size=16,fontweight='bold')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%d-%b"))
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    ax2.set(xlabel='Time Stamp', ylabel='SSH (m)')
    ax1.yaxis.set_label_coords(-0.075, 0.5)
    ax2.yaxis.set_label_coords(-0.075, 0.5)
    for item in ([ax1.xaxis.label, ax1.yaxis.label]):
        item.set_fontsize(16)
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(14)
    for item in ([ax2.xaxis.label, ax2.yaxis.label]):
        item.set_fontsize(16)
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(14)
    fig.suptitle(mdates.DateFormatter("%Y-%m-%d %H:%M")(data['time'][timeidx]),
                fontsize='20',fontweight='bold',y=0.92)
    plt.savefig('TSTiles/tile' + str(timeidx) + '.png', dpi=300)
    plt.show()
