# radar-validation
Tools to validate radar data against FVCOM (model) and ADCP (observational) data

**get_data.py**
Loads files and consolidates them into a common structure. Contains three functions:
(1) getFVCOM: Takes .p file containing model data and grabs u,v,speed, and time at each lat/lon point
(2) getRadar: Takes .p file containing radar data and grabs u,v,speed, and time at each lat/lon point
(3) getADCP: Takes .mat file containing ADCP data at one lat/lon point, must supply location in lat/lon

**do_utide.py**
Contains solveCoefs function which takes in radar or model data structure and returns the structure with a field containing u-tide coefficients at specified locations. By default, will calculate coefficients at all grid locations. Can also specify locations to calculate u-tide at (in lat,lon format) if you don't want to do all of them. 

**val_plots.py**

**calcs.py**

**error_metrics.py**

**********************************************************************************************************************************************************
**get_receiver_coefs.py**

**dictToCSV.py**

**mk_grid.py**

**mk_samples.py**


