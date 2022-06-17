# radar-validation
Tools to validate radar data against FVCOM (model) and ADCP (observational) data

**get_data.py**

Loads files and consolidates them into a common structure. Contains three functions:

    (1) getFVCOM: Takes .p file containing model data and grabs u,v,speed, and time at each lat/lon point
    
    (2) getRadar: Takes .p file containing radar data and grabs u,v,speed, and time at each lat/lon point
    
    (3) getADCP: Takes .mat file containing ADCP data at one lat/lon point, must supply location in lat/lon

**do_utide.py**

Contains solveCoefs function which takes in radar or model data structure and returns the structure with a field containing u-tide coefficients at specified locations. By default, will calculate coefficients at all grid locations. Can also specify locations to calculate u-tide at (in lat,lon format) if you don't want to do all of them. 

**calcs.py**

Does calculations of ssh gradient, vorticity, divergence, utide residuals on the radar grid scale. Contains four functions: 

    (1) vorticity: takes in a dictionary (from getFVCOM() or getRadar() functions in get_data), calcuate 
    
    (2) divergence: grab u-tide coefficients and use them to reconstruct time-series at each location
    
    (3) sshgrad: calculates spatial gradient of sea surface height (SSH) between grid points
    
    (4) residuals: from u_tide results, get residuals and calculate vorticity/divergence

**error_metrics.py**

Error metrics used for radar validation. Calculates RMSE & NRMSE of speed metric, complex correlation, and phase angles. 

**val_plots.py**
Functions to make plots with radar and model

Provides functions:

    - hex_to_rgb: Converts hex codes to rgb for plotting w pyplot 
    
    - rgb_to_dec: Converts rgb to decimal for plotting w pyplot
    
    - get_continuous_cmap: creates custom colourmap given input colours in dec
    
    - plotTS: Plots time series of input variable(s) at one location for radar OR model 
    
    - compareTS: Plots two time series: one of radar and one of model given input variable(s) 
    
    - spatialplot: Creates spatial plot of input parameter for radar or model data over whole grid 
    
    - radarcount: Creates spatial plot of number of radar measurements that pass ellipse QC test 
    
    - plotellipse: Plots values of ellipse parameter (error measure) given input radar data
    
**********************************************************************************************************************************************************
Supplemental Files (Not the Cleanest Files) 

**get_receiver_coef.py**

Loops through csv of receiver locations and dates and calculates coeffs for each

**dictToCSV.py**

Converts dictionary to csv (for CB)

**mk_grid.py**

Calculates flow metrics for all locations on radar grid for entire time period of interest (2011-2021), exports data as csv (very slow - best to run on server!)

**mk_samples.py**

Calculates flow metrics for receiver locations during deployment durations of each receiver, exports data as csv


