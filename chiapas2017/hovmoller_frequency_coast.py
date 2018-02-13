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
fout='/Users/dmelgar/Tsunamis/tehuantepec_48hrs/_output/_gauges_freq.mseed'
contour=genfromtxt('/Users/dmelgar/DEMs/SRTM15/tehuantepec400.xyz')

gmt_file='/Users/dmelgar/code/GMT/tehuantepec/hovemoller.xyz'

# gauge files
path='/Users/dmelgar/Tsunamis/tehuantepec_48hrs/_output/'
#path='/Users/dmelgar/Tsunamis/tehuantepec_mentawai/_output/'
#path='/Users/dmelgar/Tsunamis/tehuantepec_iquique/_output/'
#path='/Users/dmelgar/Tsunamis/tehuantepec_illapel/_output/'
files=glob(path+'gauge02*')

#Absolute timevector
#tout=arange(0,86300,60.)
tout=arange(0,24*3600,60.)

#Minimum amr level
min_amr=1

#Shoal depth?
H0=5.0

#Coordinates of gauges
xyz=genfromtxt('/Users/dmelgar/DEMs/SRTM15/tehuant_coast_points_filtered.xyz')

#Read all data, interpolate, put in Stream
st=Stream()
for k in range(len(files)):
    
    print '%d of %d' % (k,len(files))

    shoal=(abs(xyz[k,2])/H0)**0.25    
            
    g=genfromtxt(files[k])

    #Amr levels
    levels=g[:,0]
    
    #Time and amplitudes
    t=g[:,1]
    eta=g[:,5]*shoal
    
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

#What are the locations of the gauges?
#lon_gauge=zeros(len(files))
#lat_gauge=zeros(len(files))
#for k in range(len(files)):
#    f=open(files[k])
#    line=f.readline()
#    f.close()
#    x=float(line.split()[4])
#    y=float(line.split()[5])
#    lon_gauge[k]=x
#    lat_gauge[k]=y

lon_gauge=xyz[:,0]
lat_gauge=xyz[:,1]

#onwards, what longitudes do you want outout for
lon=arange(lon_gauge.min(),lon_gauge.max(),0.01)

#Interpolate to latitudes as well, meed this for mapproject in gmt
f=interp1d(lon_gauge, lat_gauge, kind='linear')
lat=f(lon)

#how many times do we have?
t=fobs.copy()

#output for havemoller
eta=zeros((len(lon),len(t)))

for ktime in range(len(t)):
    
    #form vecotr of etas at this time slice
    eta_slice=zeros(len(lon_gauge))
    for kgauge in range(len(lon_gauge)):
        eta_slice[kgauge]=st[kgauge].data[ktime]
    
    #Interpolate
    f=interp1d(lon_gauge, eta_slice, kind='linear')
    eta_interp=f(lon)
    
    #Place in output varible
    eta[:,ktime]=eta_interp
    

#Get distance from coast point to edge of shelf "shelf width")
p=Geod(ellps='WGS84')
f=interp1d(lon_gauge, lat_gauge, kind='linear')
lat=f(lon)
dist_out=zeros(len(lon))
for k in range(len(lon)):
    az,baz,dist=p.inv(lon[k]*ones(len(contour)),lat[k]*ones(len(contour)),contour[:,0],contour[:,1])
    dist=dist/1000
    dist_out[k]=dist.min()
    
#Write to  \a file
out=c_[lon,lat,dist_out]
#savetxt(u'/Users/dmelgar/code/GMT/tehuantepec/shelf.txt',out,fmt='%.4f')


#plotaroo    
fig = plt.figure(figsize=(7,7))
gs = gridspec.GridSpec(2, 1,height_ratios=[1,3])
ax2=plt.subplot(gs[0])
ax1=plt.subplot(gs[01])


#Meshgrid tsunami
Y,X=meshgrid(t*3600.,lon)
lon_out=X.ravel()
time_out=Y.ravel()
eta_out=eta.ravel()

#Meshgrid lon,lat
Y,X=meshgrid(t*3600.,lat)
lat_out=X.ravel()

ax1.set_yscale('log')
o=ax1.pcolormesh(lon,t[1:]*3600,log10(eta[:,1:].T),cmap=plt.cm.magma,vmin=1)
ax1.set_xlabel('Longitude')
ax1.set_ylabel('CPH')
ax1.set_ylim([1e-1,1e1])
#cb=fig.colorbar(o)
#cb.set_label('PSD')

#gauges
ax1.plot([-96.1169,-96.1169],[1e-1,1e1],'--',c='b')
ax1.plot([-92.4012,-92.4012],[1e-1,1e1],'--',c='b')
ax1.plot([-96.4917,-96.4917],[1e-1,1e1],'--',c='b')
ax1.plot([-95.1968,-95.1968],[1e-1,1e1],'--',c='b')

ax2.plot(lon,dist_out/10,lw=2,c='k')
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Width (km)')

plt.show()
#savetxt(gmt_file,c_[lon_out,lat_out,time_out,eta_out],fmt='%.4f')
