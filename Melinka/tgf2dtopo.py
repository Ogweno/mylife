from mudpy import inverse
from obspy.core import UTCDateTime


project_name='Melinka_usgs'
model_name=u'maule.mod'
tgf_file='seafloor.sta'
coast_file=None
fault_name=u'melinka_usgs.fault'
time_epi=UTCDateTime('2016-12-25T14:22:26.00')
topo_dx_file=u'/Users/dmelgar/DEMs/SRTM30/melinka.dx.grd'
topo_dy_file=u'/Usersr/dmelgar/DEMs/SRTM30/melinka.dy.grd'
home='/Users/dmelgar/Slip_inv/'
tsun_dt=2.
maxt=120.

#inverse.make_tgf_dtopo(home,project_name,model_name,tgf_file,coast_file,fault_name,
#            time_epi,tsun_dt,maxt,topo_dx_file,topo_dy_file,
#            topo_effect=False,instantaneous=True,hot_start=0)
            
#inverse.tsunami_gf(home,project_name,model_name,fault_name,197)

tlims=[0.,7200.]
dt=60.
inverse.tsunami2sac(home,project_name,model_name,fault_name,tlims,dt,time_epi,0)