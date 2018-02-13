from numpy.random import uniform
from mtspec import mtspec
from obspy import read
from numpy import r_,mean,where,size,zeros,c_,ones,savetxt,arange
from netCDF4 import Dataset
from matplotlib import pyplot as plt
from string import rjust
#####Run time parameters

duration=24*3600 #in seconds
dt=60*15
fcorner=[1./(2*3600)]   #defines the tsunami band
month=['06','06','06','07']
darts=['46404','46407','46411','46419']
dart_depths=['2736','3268','4286','2773']
path_darts='/Users/dmelgar/tidegauge_noise/DART/sac/'
max_eta=0.004
bathy_file=u'/Users/dmelgar/DEMs/crescent_city/crescent15_big.grd'
fout='/Users/dmelgar/tidegauge_noise/data/cres/_dtopo/time_varying_24hrs.dtopo'
#####


Ntimeseries=int(duration/dt)


def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt




#read them all
for k in range(len(darts)):
    sta=darts[k]
    st=read(path_darts+sta+'_2017_'+month[k]+'.sac')
    if k==0:
        data=st[0].data
        data=data-mean(data)
    else:
        data=r_[data,st[0].data-mean(st[0].data)]




data=highpass(data,fcorner,1./(15*60),2)
i=where(abs(data)>max_eta)[0]
data[i]=max_eta

ocean_noise=data.copy()
Nsamples=len(ocean_noise)



#read bathymetry file
grd = Dataset(bathy_file, 'r', format='NETCDF4')
x=grd.variables['lon'][:]
y=grd.variables['lat'][:]
z=grd.variables['z'][:]


#write it
decimate=20
decimate_x=range(0,len(x),decimate)
decimate_y=range(0,len(y),decimate)
coords=zeros((len(decimate_x)*len(decimate_y),2))
zref=zeros(len(decimate_x)*len(decimate_y))
row=0
y=y[::-1]

for klat in decimate_y:
    for klon in decimate_x: 
        coords[row,:]=r_[x[klon],y[klat]]
        zref[row]=z[len(y)-klat-1,klon]
        row+=1


zout=zeros((len(zref),Ntimeseries))
indices=uniform(0,Nsamples-Ntimeseries-1,size(zref))
indices=indices.astype('int')

for kz in range(len(zref)):
    zout[kz,:]=ocean_noise[indices[kz]:indices[kz]+Ntimeseries]


    
#find land points
land=where(zref>0)[0]   
             
out=c_[zeros(len(zref)),coords,zeros(len(zref))]
for kz in range(Ntimeseries):
    time=dt*(kz+1)
    eta=zout[:,kz]
    eta[land]=0
    out_slice=c_[time*ones(len(zref)),coords,eta]
    out=r_[out,out_slice]
    
    

    
##plot it?
#for k in range(Ntimeseries-1):
#    print k
#    plt.figure()
#    i=arange(len(coords))+k*len(coords)
#    s=plt.scatter(out[i,1],out[i,2],c=out[i,3]*1000,lw=0,cmap=plt.cm.seismic,vmin=-4,vmax=4)
#    cb=plt.colorbar(s)
#    cb.set_label('Sea surface height (mm)')
#    plt.title('t = %.2fhrs' % (k*dt/3600.))
#    num=rjust(str(k),4,'0')
#    plt.savefig('/Users/dmelgar/tidegauge_noise/plots/dtopo/'+num+'.png')
#    plt.close()
#    


savetxt(fout,out,fmt='%.1f\t%.8f\t%.8f\t%.6f')








