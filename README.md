# OceanDiv v0.1 

OceanDiv is a tool for quantifying smoothly varying ocean heat transport divergences (convergences) from estimates of ocean heat content, net surface heat fluxes, and their associated uncertianties. 

## Methodology
OceanDiv uses a forward Kalman filter and backward Rauch–Tung–Striebel smoother to estimate ocean heat transport convergences (HTC) from observations of full-column ocean heat content (OHC) and net surface heat fluxes. A [Kaman filter](https://en.wikipedia.org/wiki/Kalman_filter) is a way of combining predictions from a model with uncertain observations to generate an improved state estimate with associated uncertainties. The use of the RTS smoother allows each estimate of HTC(t) using all available (past and future) OHC and surface flux observations. This methodology is similar to that employed by Kelly et al. ([2014](http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-12-00131.1)) with some key differences in the treatment of uncertainties and the way that HTC is explitly incorporated into the Kalman filter state transition matrix. In order to run OceanDiv, it is necessary to provide at least two different estimates of each observed quantity (ocean heat content and net surface heat flux) as the discrepancies between different products are used to estimate observational uncertainties. 


## User Guide
OceanDiv is invoked by running the run_oceandiv.py script with a path towards an appropriately configured namelist file. 

``` > python2.7 run_oceandiv.py config/namelist.ini ```

### Dependencies
You will need a python2.7 enviroment with the following packages (and their associated dependencies) installed: ``` iris, numpy, ConfigParser, argparse, and scipy```.


### Configuring the namelist
The namelist file contains paths to datasets and various options required for OceanDiv to run and is configured into sections using the [ini format] (https://en.wikipedia.org/wiki/INI_file). The annotated example below shows how the namelist should be configured. 

```
[metadata]
# IMPORTANT - all datasets must have the same grid coordinates for OceanDiv to work. 
nohc = 2 # Number of different ocean heat content data sets (must be > 1).
nflx = 2 # Number of net surfacd heat flux data sets (must be > 1). 
nx = 360 # Number of grid points in x-direction. 
ny = 180 # Number of grid points in y-direction.
nt = 336 # Number of data points in time axis.  

[output]
dir = /Users/chris_roberts/data/energy_flows/kalman_output/ # Location to save output.
name = global_1x1  # Experiment name to use in output files.
print_stdout = True # Boolean flag used to switch on/off output to the terminal. 

[basins]
dir = /Users/chris_roberts/data/energy_flows/ # Directory containing basin masks. 
f = 1x1_grid_mask_r360x180.nc # 2D netcdf file with each cell containing a value corresponding to a specifed ocean basin. 
var = grid_mask # Name of the 2D netcdf variable containing the basin masks. 
mdi = -1 # Missing data indicator used for the basin masks.  

[areas]
dir = /Users/chris_roberts/data/energy_flows/ # Directory containing grid area weights.
f = cell_area_r360x180.nc # Netcdf file containing grid area weights.
var = cell_area # Name of variable
mdi = -1 # Missing data indicator. 
calc_basin_totals = False # Boolean flag to switch on/off use of area totals. Default is to use area means. 

[kfilter]
smooth = True # Boolean flag that switches on/off application of a low-pass Butterworth filter prior to Kalman filter.
cutoff = 12 # Time window used for low-pass filter. 
dt = 2629800.0 # Duration (in seconds) of each data point (2629800.0  = 1 month).  
oht_error_scaling = 2 # Scale factor used to estimate uncertainty in persistence of heat
                      # transport convergence as a multiple of the uncertainty in the persistence
                      # of net surface fluxes (default = 2). 

[ohc1] # Details of OHC data set number 1
dir = /Users/chris_roberts/data/energy_flows/ # Path to data
f = EN.4.1.1.f.analysis.g10.ocean_heat_content_0-10000m_mon_deseasoned_vs_2006-2012_r360x180.nc # NetCDF file - units of J/m2
var = ocean_heat_content # variable name
t0 = 120 # Inital time index. 
mdi = 1e20 # Missing data indicator

[ohc2] # Details of OHC data set number 2, same as above. 
dir = /Users/chris_roberts/data/energy_flows/
f = pot_anal_ocean_heat_content_0-10000m_mon_deseasoned_vs_2006-2012_r360x180.nc
var = ocean_heat_content
t0 = 120
mdi = 1e20

[flx1] # Details of flux data set number 1. 
dir = /Users/chris_roberts/data/energy_flows/update_2016-02-18/  # Path to data
f = fmass_eraint_mon_deseasoned_vs_2006-2012_r360x180_spherical.nc # NetCDF file - units of W/m2
var = Fmass # variable name
t0 = 0 # Initial time index
mdi = 1e20  # Missing data indicator. 

[flx2] # Details of flux data set number 1, same as above.
dir = /Users/chris_roberts/data/energy_flows/update_2016-02-08/
f = fmass_merra_mon_deseasoned_vs_2006-2012_r360x180_spherical_masked.nc
var = Fmass
t0 = 0
mdi = 1e20




