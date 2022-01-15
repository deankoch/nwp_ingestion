# Install pkern from Python
# Dean Koch
# January 11, 2022
# 
#
# INTRODUCTION
#
# This script documents how I installed and set up an R
# environment from Python, using rpy2. The purpose is to
# access the pkern R package for simple kriging from Python.
# 
# I've switched from Windows 10 to Linux Mint v20.3 to work on
# this code, since rpy2 doesn't have good windows support at the
# moment. As with the rest of the NWP_ingestion code, I'm using
# the "herbie" anaconda environment defined in /environment.yml.
#
# rpy2 is an interface to R functions and objects from Python. It
# comes with its own portable build of R, so there is no need
# to worry about installing R itself on your system, but we do
# need to install a few R packages in order to access my kriging
# code. 
#
#
# INSTALL ME FIRST
#
# Some of the R packages required below have their own external
# dependencies. They won't install unless the following system
# libraries are accessible from rpy2's embedded R enviroment:
#  
# libcurl, udunits2, libgdal, geos, and proj.4
#
# If you run into R package installation errors related to these
# dependencies, you may need to tinker with R's environmental
# variables $PATH and $PKG_CONFIG_PATH.
# 
# In my case R couldn't find libcurl and udunits2 until I
# installed them into my conda environment and pointed R to
# the /lib folder for that environment (see below).
# 
#
# INSTALL RPY2
#
# rpy2 for Python 3.9 can be installed by activating your conda
# environment in a terminal and running:
# 
#   conda install -c conda-forge rpy2
#
# It was necessary to switch to the "conda-forge" channel (instead
# of "r" channel) to get a build that works with Python 3.9
#
#
# INITIALIZE R FROM PYTHON
# 
# R is initialized when we import the "robjects" module in Python:

# rpy2 modules
import rpy2.robjects as R
import rpy2.robjects.packages as Rlib


# %%
#
'''------- Set up paths --------'''

# The following directories should be replaced by the user's own:

# path to R helper function definitions
nwp_dir = '/mnt/e/coding-projects/nwp_ingestion'
rhelp_path = nwp_dir + '/pkern_with_rpy2/R_helpers.R'

# path to /lib for the anaconda environment in use  
lib_path = '/mnt/e/anaconda3/envs/herbie/lib'

#
# %%
#
'''------- Set up the R environment --------'''

# We will install the required R packages from inside R, using
# a few helper R functions defined in the "R_helpers.R" file.
# The way this works is we use Python's open(...).read() to load
# the R script text and pass it to R via R.r(...)

# load R helper functions in embedded R
R.r(open(rhelp_path).read())
# these R functions callable from Python now! eg "set_env" below

# add conda library location to $PATH and $PKG_CONFIG_PATH
R.r['set_ev'](lib_path, lib_path + '/pkgconfig')


# %%
'''------- install R dependencies --------'''

# This last chunk installs the following R packages:
# 'raster', 'sf', 'codetools', 'pkern' (and their dependencies)

# import R base functions and set CRAN mirror for downloads
Rutils = Rlib.importr('utils')
Rutils.chooseCRANmirror(ind=0)

# The calls below will download and compile from source a number
# of packages. Some of these may take some time to complete

# install spatial libraries for R and a suggested package for rpy2
pkg_tocheck = ('raster', 'sf', 'codetools')
for pkg in pkg_tocheck:
    if not Rlib.isinstalled(pkg):
        Rutils.install_packages(pkg)

# pkern is not yet on CRAN, but we can fetch it directly from
# github using an additional library, devtools

# install, then import devtools (remove when pkern gets on CRAN)
if not Rlib.isinstalled('devtools'):
    Rutils.install_packages('devtools')
Rdevtools = Rlib.importr('devtools')

# install pkern via devtools
if not Rlib.isinstalled('pkern'):
    Rdevtools.install_github('deankoch/pkern')

# Once a package is installed it will be accessible in
# subsequent R sessions, so you should only have to run this
# script once.
#
# if any of these installations fail, fix the problem (system
# libraries not on R's $PATH?) then re-run this script.
# install_packages should detect when a package is already
# installed and skip it.