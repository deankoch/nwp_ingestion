# Simple downscaling with pkern
# Dean Koch
# January 14, 2022
#
# OVERVIEW
#
# This script interpolates temperature data onto a higher-resolution
# target grid by running a pkern-based R workflow from Python.
#
# The input files target_grid_simple.nc and temperature_series.nc should
# be found in testdata/pkern_examples (run the previous two scripts
# in this series to make them).
#
# The script writes netcdf output to temperature_series_output.nc in
# testdata/pkern_examples, along along with several png figures in the
# new folder pkern_graphics

import os
import numpy as np
import geopandas
import xarray as xr
import matplotlib.pyplot as plt

# rpy2 modules
import rpy2.robjects as R
import rpy2.robjects.packages as Rlib

# %%
'''------- define working directory and filepaths --------'''

# path to data subdirectory for output
nwp_dir = '/mnt/e/coding-projects/nwp_ingestion'
pkern_data_dir = nwp_dir + '/testdata/pkern_examples'

# path to R helper function definitions
rhelp_path = nwp_dir + '/pkern_with_rpy2/R_helpers.R'

# paths to input data, target grid, output data
krig_input_path = pkern_data_dir + '/temperature_series.nc'
target_grid_path = pkern_data_dir + '/target_grid_simple.nc'
krig_output_path = pkern_data_dir + '/temperature_series_output.nc'

# new graphics directory
pkern_graphics_dir = nwp_dir + '/pkern_graphics' 
if not os.path.isdir(pkern_graphics_dir):
    print('creating data directory at: ' + pkern_graphics_dir)
    os.mkdir(pkern_graphics_dir)

# paths to output graphics
img_vario_path = pkern_graphics_dir + '/variogram_simple.png'
img_krig_result_path = pkern_graphics_dir + '/kriging_result_simple.png'

# %%
'''------- Set up R to run pkern --------'''
#
# R was initialized by the earlier rpy2 import call. To proceed we
# need to load some required libraries and helper functions
 
# load R helper functions
R.r(open(rhelp_path).read())

# load the "pkern", "sf", "raster" R packages 
Rsf = Rlib.importr('sf')
Rraster = Rlib.importr('raster')
Rpkern = Rlib.importr('pkern')

# %%
'''------- convert target grid to R-friendly objects --------'''
#
# R integration with rpy2 is seamless when it comes to simple objects
# like vectors. More complicated objects like a rasters and points
# datasets aren't directly readable in R from python so we have to first
# extract their essential components as strings and vectors.
#
# In this case I've defined the target grid via a netcdf file, but if we
# want it could also be constructed directly from input parameters like
# dimension, extent, etc. The chunk below extracts the relevant grid info: 

# open the netcdf file
target_grid = xr.open_dataset(target_grid_path)

# extract projection, grid line locations, and target resolution
target_crs = target_grid['spatial_ref'].crs_wkt
target_glx = target_grid['x'].to_numpy()
target_gly = target_grid['y'].to_numpy()
target_res = [target_gly[1] - target_gly[0], target_glx[1] - target_glx[0]]

# in R: define a named list containing grid line locations
Rgyx = R.r['list'](R.FloatVector(target_gly), R.FloatVector(target_glx))
Rgyx = R.r['setNames'](Rgyx, nm=R.StrVector(['y', 'x']))

# in R: append grid dimensions and set element names expected by pkern
Rgdim = R.r['sapply'](Rgyx, 'length')
Rgnames = R.StrVector(['gdim', 'gyx', 'gres'])
Rg = R.r['list'](Rgdim, Rgyx, R.FloatVector(target_res))
Rg = R.r['setNames'](Rg, nm=Rgnames)
# "Rg" is (the python interface to) an R list defining the target grid
  

# %%
'''------- convert input points to R-friendly objects --------'''
# 
# In this example, the input data happens to be arranged in a grid;
# However this is not a requirement the pkern algorithm. All we need
# are a set of coordinates, their projection info, and the associated
# data values. 
#
# The following chunk reshapes the input dataset, making a points
# geometry object (geopandas) along with a list of data vectors
# (one per time slice).

# open the temperature data
input_data = xr.open_dataset(krig_input_path)
input_crs = input_data['spatial_ref'].crs_wkt
input_nx = len(input_data['x'])
input_ny = len(input_data['y'])
input_ntime = len(input_data['time'])

# expand grid line positions to get coordinates of all grid points
input_x = input_data['x'].to_numpy().repeat(input_ny)
input_y = np.tile(input_data['y'].to_numpy(), input_nx)
input_time = input_data['time']

# create a points geometry object defining temperature data locations
pts = geopandas.points_from_xy(input_x, input_y, crs=input_crs)

# Note that it's important that the projection of "pts" (above)
# matches the projection of "target_grid". If they have different
# projections, you must transform the input data to the same projection
# (crs) as the target grid!

# transform point locations to target projection
pts_t = pts.to_crs(target_crs)

# In this example the transformation does nothing, since input_crs and
# target_crs are identical. However in typical usage the target grid
# and input points will have different projections, so the above step
# produces new coordinate values for pts:

# in R: construct a dataframe and CRS object from transformed coords
Rcrs = R.r['st_crs'](target_crs)
Rpts = R.r['data.frame'](R.FloatVector(pts_t.x), R.FloatVector(pts_t.y))
Rpts = R.r['setNames'](Rpts, nm=R.StrVector(['x', 'y']))

# in R: construct point geometry locating the observed data
Rsf = R.r['st_as_sf'](Rpts, coords=R.StrVector(['x', 'y']), crs=Rcrs)


# %%
'''------- convert input data to R-friendly object --------'''
#
# The data (temperature) values themselves are unaffected by projections.
# But we need to reshape them into a simpler form to pass to R. This
# chunk vectorizes each time slice, building a list of data vectors in R. 

# unpack data to list of vectors (in column-vectorized order)
vname = list(input_data.data_vars)[0]
pts_data = list(input_data[vname].to_numpy())
pts_data_vectors = [x.flatten(order='F') for x in pts_data]

# in R: make list of data vectors, one per time slice
keys = range(len(pts_data_vectors))
pts_data_list = [R.FloatVector(pts_data_vectors[x]) for x in keys]
pts_data_dict = {'d' + str(x) : pts_data_list[x] for x in keys}
Rpts_data = R.ListVector(pts_data_dict)


# %%
'''------- snap input points to grid --------'''
#
# All the input data are now loaded in R. We begin the workflow
# by snapping the input points to the target grid

# in R: snap points to regular grid without any duplication
Rgsnap = R.r['pkern_snap'](Rg, pts=Rsf)
Rsgdim = Rgsnap.rx2('sg').rx2('gdim')

# in R: data can now be mapped to the grid using a helper function
Rsgdata = R.r['get_sglist'](Rpts_data, Rgsnap)


# %%
'''------- fit the variogram --------'''

# the variogram is calculated over all (time) layers

# in R: compute sample variograms
Rvario = R.r['pkern_vario'](Rgsnap, Rsgdata)

# uncomment this line to plot the sample variogram point cloud
# R.r['pkern_vario_plot'](Rvario, Rpvario)

# At small distances these look very similar to theory. However
# in many cases the variogram shows a big upward swing at the tails,
# whereas theory predicts a flat approach to the sill value. Since
# we can't model the unexpected tail behaviour, we exclude the tails
# from model fitting by setting a maximum distance "dmax" in the next
# function call (with units meters).   

# in R: fit a theoretical Matern X Matern variogram, excluding lags > 25km
dmax = 25 * 1000
Rpvario = R.r['pkern_vario_fit'](
    Rvario, 
    ypars='mat',
    xpars='mat', 
    dmax=dmax) 

# report fitted parameters then (in R) print variogram results to png file
print(Rpvario)
plotpars = R.r['setNames'](R.r['list'](dmax), 'dmax')
R.r['png'](img_vario_path, height=600, width=800, units='px')
R.r['pkern_vario_plot'](Rvario, Rpvario, plotpars=plotpars)
R.r['dev.off']()


# %%
'''------- run the kriging algorithm --------'''

# in R: krig on all time slices, convert output directly to numpy array
krig_result = np.array(R.r['run_krig'](Rsgdata, Rgsnap, Rpvario))

# reshape into vectorization order expected by xarray
zpred = []
for res in krig_result:
    zpred.append(res.reshape(Rgdim, order='F'))

# build xarray dataset from the output
output_data = xr.Dataset(
    data_vars=dict(temperature=(['time', 'y', 'x'], np.asarray(zpred))),
    coords=dict(
        x=target_glx,
        y=np.flip(target_gly),
        time=input_time.to_numpy())
    )

# write the projection info
output_data = output_data.rio.write_crs(target_crs)

# save output as netcdf
output_data.to_netcdf(path=krig_output_path)
#
# %%
'''------- plot the result for first time slice --------'''

# initialize a two-pane plot
f, (ax1, ax2) = plt.subplots(1, 2)

# plot the source data on the left
input_data[vname][0,:,:].plot(ax=ax1)
ax1.set_title('source data')
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)

# plot kriging output on the right
output_data['temperature'][0,:,:].plot(ax=ax2)
ax2.set_title('kriging output')
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

# save to png
plt.tight_layout()
plt.savefig(img_krig_result_path, facecolor='white', dpi=400)


# %%
