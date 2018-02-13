from numpy import zeros,genfromtxt,arange
from mudpy import inverse

from mudpy import runslip
import numpy as np
from obspy.core import UTCDateTime
########                            GLOBALS                             ########
home='/Users/dmelgar/Slip_inv/'
project_name='Nepal_ALOS'
run_name='ALOS_LRGPS_gpsonly'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####close   

init=0 #Initalize project
make_green=0 #Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
G_from_file=0 # =0 read GFs and create a new G, =1 load G from file
invert=1  # =1 runs inversion, =0 does nothing
###############################################################################

###############          view  Green function parameters               #############
ncpus=1
hot_start=0  #Start at a certain subfault number
model_name='avouac.mod'   #Velocity model
fault_name='nepal_10.fault'    #Fault geometry
GF_list='nepal_lrgps_t048.gflist'#What GFs are to be computed for each station
#############               Inversion Parameters               ################# 
time_epi=UTCDateTime('2015-04-25T06:11:26')
epicenter=np.array([84.708,28.147,15]) 
rupture_speed=2.7 #Fastest rupture allowed in km/s
num_windows=1
reg_spatial=np.logspace(-6,1,num=30) #Set to False if you don't want to use it
reg_temporal=None#np.logspace(-6,-6,num=1) #Set to False if don't want to use it
nstrike=20 ; ndip=15 ; nfaults=(nstrike,ndip) #set nstrike to total no. of faults and ndip to 1 if using Tikh
beta=45 #Rotational offset (in degrees) applied to rake (0 for normal)
Ltype=0 # 0 for Tikhonov and 2 for Laplacian
solver='nnls' # 'lstsq','nnls'
top='free' ; bottom='locked' ; left='locked' ; right='locked' #'locked' or 'free'
bounds=(top,bottom,left,right)
################################################################################e=
########      Run-time modifications to the time series             ############
weight=False
decimate=None  #Decimate by constant (=None for NO decimation)
bandpass=None #Corner frequencies in Hz =None if no filter is desired
################################################################################
G_name='test'


m=zeros((600,1))
iss=arange(0,len(m),2)
ids=arange(1,len(m),2)
mod=genfromtxt(u'/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/output/inverse_models/models/ALOS_GPS_3.3_20win_pulse_vall3.0001.inv.total')
m[iss,0]=mod[:,8]
m[ids,0]=mod[:,9]
mrot=inverse.ds2rot(m,45)
G=inverse.getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
            rupture_speed,num_windows,decimate,bandpass)
d=inverse.getdata(home,project_name,GF_list,decimate,bandpass=None)
ds=G.dot(mrot)

