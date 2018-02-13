from mudpy import inverse
from obspy import UTCDateTime
from numpy import array

home='/Users/dmelgar/Slip_inv/'
project_name='Chiapas_hernandez_new'
run_name='tsun'
model_name='hernandez.mod'
tgf_file='seafloor.gflist'
coast_file=None
time_epi=UTCDateTime('2017-09-08T04:49:17.3')
epicenter=array([-94.103,14.761,45.9])
#parameters of the dispalcement synthetics ALREADY made
tsun_dt=5.0
maxt=320
topo_dx_file='/Users/dmelgar/DEMs/SRTM30/tehuantepec2m.grd'
topo_dy_file='/Users/dmelgar/DEMs/SRTM30/tehuantepec2m.grd'
fault_name='chiapas_new.fault'
tlims=[0,3600*2.0]
dt=60.0
hot_start=0


#This makes the dtopo

#inverse.make_tgf_dtopo(home,project_name,model_name,tgf_file,coast_file,fault_name,
#            time_epi,tsun_dt,maxt,topo_dx_file,topo_dy_file,
#            topo_effect=False,instantaneous=True,hot_start=0,average_instantaneous=4)
            
            
#This runs geoclaw
#hot_start=0
#inverse.tsunami_gf(home,project_name,model_name,fault_name,hot_start)       

###Make SAC files     
inverse.tsunami2sac(home,project_name,model_name,fault_name,tlims,dt,time_epi,hot_start)     
