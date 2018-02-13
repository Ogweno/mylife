from numpy import genfromtxt,where,r_,savetxt

dtopo_file= '/Users/dmelgar/Maule2010/tsunami/dtopo/gps_tg_kinematic.dtopo'
dtopo_file_static= u'/Users/dmelgar/Maule2010/tsunami/dtopo/gps_tg_static.dtopo'

dtopo=genfromtxt(dtopo_file)

i0=where(dtopo[:,0]==0)[0]
tmax=dtopo[:,0].max()
imax=where(dtopo[:,0]==tmax)[0]
dtopo[imax,0]=1
static_dtopo=r_[dtopo[i0,:],dtopo[imax,:]]
savetxt(dtopo_file_static,static_dtopo,fmt='%d\t%13.8f\t%13.8f\t%.6e')