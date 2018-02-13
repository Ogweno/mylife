from obspy import read
from scipy.interpolate import interp1d
from mudpy.forward import lowpass
from glob import glob
from numpy import arange

files=glob('/Users/dmelgar/Slip_inv/Chiapas_hernandez/data/waveforms/*.disp*')
dt=0.25
fcorner=0.49

for k in range(len(files)):
    
    print files[k]
    
    st=read(files[k])
    dt_data=st[0].stats.delta
    tmax=dt_data*st[0].stats.npts-dt_data
    t=arange(0,tmax,dt)
    
    #filter
    d=lowpass(st[0].data,fcorner,1./dt_data,2,zerophase=True)
    
    #Interpolate
    f=interp1d(st[0].times(),d)
    dinterp=f(t)
    
    #filter again to be safe about upsampled data
    dout=lowpass(dinterp,fcorner,1./dt,2,zerophase=True)
    
    #metadata
    st[0].stats.delta=dt
    
    #output
    st[0].data=dout
    st.write(files[k],format='SAC')