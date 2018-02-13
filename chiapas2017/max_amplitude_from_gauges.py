from numpy import genfromtxt,where,savetxt,zeros
from glob import glob

path=u'/Users/dmelgar/Tsunamis/tehuantepec_hf_all/_output/'
files=glob(path+'gauge9*')

#xyz=genfromtxt('/Users/dmelgar/DEMs/SRTM15/tehuant_coast_points_filtered.xyz')
min_amr=3
max=3.5

#Shoal depth?
H0=2.0
xyz=genfromtxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_domain_points_inwater.xyz')
out=xyz.copy()
for k in range(len(files)):
    
    print files[k]
    
    #read gauge
    g=genfromtxt(files[k],usecols=[0,5])
    i=where(g[:,0]>=min_amr)[0]
    
    #get shoaling correction
    #shoal=(abs(xyz[k,2])/H0)**0.25
    shoal=1
    out[k,2]=g[i,1].max()*shoal
    
i=where(out[:,2]<max)[0]
out=out[i,:]
    
savetxt(u'/Users/dmelgar/code/GMT/tehuantepec/max_eta_all.xyz',out,fmt='%.4f')