from numpy.random import uniform
from mtspec import mtspec
from obspy import read
from numpy import r_,mean,where,size,zeros,c_,ones,savetxt
from netCDF4 import Dataset
#####Run time parameters

duration=48*3600 #in seconds
dt=60.
fcorner=[1./(2*3600)]   #defines the tsunami band
month=['06','06','06','07']
darts=['46404','46407','46411','46419']
dart_depths=['2736','3268','4286','2773']
path_darts='/Users/dmelgar/tidegauge_noise/DART/sac/'
max_eta=0.004
bathy_file=u'/Users/dmelgar/DEMs/crescent_city/crescent15_big.grd'
fout='/Users/dmelgar/tidegauge_noise/data/cres/_dtopo/noise_from_DARTs.dtopo'
#####





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

indices=uniform(0,Nsamples,size(z))
indices=indices.astype('int')

zout=ocean_noise[indices]

#if on alnd don;t move it
land=where(z.ravel()>0)[0]
zout[land]=0

#write it
decimate=20
decimate_x=range(0,len(x),decimate)
decimate_y=range(0,len(y),decimate)
coords=zeros((len(decimate_x)*len(decimate_y),2))
row=0
y=y[::-1]

for klat in decimate_y:
    for klon in decimate_x: 
        coords[row,:]=r_[x[klon],y[klat]]
        row+=1
zout=zout[0:len(coords)]
       
out_time=r_[zeros(len(zout)),ones(len(zout))]
out_eta=r_[zeros(len(zout)),zout]        
out_coords=r_[coords,coords]
out=c_[out_time,out_coords,out_eta]

savetxt(fout,out,fmt='%.1f\t%.8f\t%.8f\t%.6f')








