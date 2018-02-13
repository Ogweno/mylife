from matplotlib import pyplot as plt
from numpy import genfromtxt,where,argmin,ones,array,r_,arange,zeros,c_,linspace,meshgrid,savetxt,log10
from pyproj import Geod
from glob import glob
from scipy.interpolate import interp1d
from string import rjust
from obspy import Stream,Trace
from mudpy import gmttools
from pyproj import Geod
from mtspec import mtspec
import matplotlib.gridspec as gridspec

cmtsun=gmttools.gmtColormap(u'/Users/dmelgar/code/python/cpt/dkbluered_modif.cpt')
peak_tsunami=12

path_out='/Users/dmelgar/Chiapas2017/plots/tsunam_coast_points/'
fout='/Users/dmelgar/Tsunamis/tehuantepec_48hrs/_output/_alldomain_freq.mseed'
contour=genfromtxt('/Users/dmelgar/DEMs/SRTM15/tehuantepec400.xyz')

gmt_file='/Users/dmelgar/code/GMT/tehuantepec/hovemoller.xyz'

# gauge files
path='/Users/dmelgar/Tsunamis/tehuantepec_48hrs/_output/'
#path='/Users/dmelgar/Tsunamis/tehuantepec_mentawai/_output/'
#path='/Users/dmelgar/Tsunamis/tehuantepec_iquique/_output/'
#path='/Users/dmelgar/Tsunamis/tehuantepec_illapel/_output/'
files=glob(path+'gauge9*')

#Absolute timevector
#tout=arange(0,86300,60.)
tout=arange(0,24*3600,60.)

#Minimum amr level
min_amr=1

#Shoal depth?
H0=5.0

#Coordinates of gauges
xyz=genfromtxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_domain_points_inwater.xyz')

#Read all data, interpolate, put in Stream
st=Stream()
for k in range(len(files)):
    
    print '%d of %d' % (k,len(files))
  
            
    g=genfromtxt(files[k])

    #Amr levels
    levels=g[:,0]
    
    #Time and amplitudes
    t=g[:,1]
    eta=g[:,5]
    
    i=where(g[:,0]>=min_amr)[0]
    t=t[i]
    g=g[i]
    
    #Keep only stuff below peak tsunamo
    #i=where(abs(eta)<=peak_tsunami)[0]
    #t=t[i]
    #eta=eta[i]
       
    #Keep only stuff above AMR 3
    #i=where(levels>=3)[0]
    #t=t[i]g
    #eta=eta[i]
    
    #resample to correct times
    f=interp1d(t, eta, kind='linear')
    eta_interp=f(tout)
    
    
    #Get spectrum
    delta=60.
    Sobs, fobs = mtspec(
        data=eta_interp, delta=delta, time_bandwidth=3.5,
        number_of_tapers=5, nfft=len(eta_interp))
        
    tr=Trace()
    tr.data=Sobs
    tr.stats.delta=delta
    st+=tr
    
st.write(fout,format='MSEED')


lon_gauge=xyz[:,0]
lat_gauge=xyz[:,1]


#Loop and look for frequencies of itnerest

f=[0.35,0.8,1.28,1.76,2.07]
f1=zeros(len(st))
f2=zeros(len(st))
f3=zeros(len(st))
f4=zeros(len(st))
f5=zeros(len(st))

for k in range(len(st)):
    
    print k
    
    i=argmin(abs(fobs*3600-f[0]))
    f1[k]=st[k].data[i]
    
    i=argmin(abs(fobs*3600-f[1]))
    f2[k]=st[k].data[i]
    
    i=argmin(abs(fobs*3600-f[2]))
    f3[k]=st[k].data[i]
    
    i=argmin(abs(fobs*3600-f[3]))
    f4[k]=st[k].data[i]
    
    i=argmin(abs(fobs*3600-f[4]))
    f5[k]=st[k].data[i]
    

savetxt(u'/Users/dmelgar/code/GMT/tehuantepec/psd_f1.xyz',c_[lon_gauge,lat_gauge,f1],fmt='%.4f')
savetxt(u'/Users/dmelgar/code/GMT/tehuantepec/psd_f2.xyz',c_[lon_gauge,lat_gauge,f2],fmt='%.4f')
savetxt(u'/Users/dmelgar/code/GMT/tehuantepec/psd_f3.xyz',c_[lon_gauge,lat_gauge,f3],fmt='%.4f')
savetxt(u'/Users/dmelgar/code/GMT/tehuantepec/psd_f4.xyz',c_[lon_gauge,lat_gauge,f4],fmt='%.4f')
savetxt(u'/Users/dmelgar/code/GMT/tehuantepec/psd_f5.xyz',c_[lon_gauge,lat_gauge,f5],fmt='%.4f')
