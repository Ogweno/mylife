from obspy import read
import glob
from datetime import timedelta
from numpy import zeros,r_
import os

f1=glob.glob(u'/Users/dmelgar/Cascadia_M9/planar_mseed/*.HNE.mseed')
f2=glob.glob(u'/Users/dmelgar/Cascadia_M9/planar_mseed/*.HNZ.mseed')
f=f1+f2
out='/Users/dmelgar/Cascadia_M9/tmp/'

noise=zeros(60*50)
for k in range(len(f)):
    
    st=read(f[k])
    #st[0].data=r_[noise,st[0].data]
    st[0].data = (st[0].data).astype('i4')
    fname=os.path.split(f[k])[1]
    st.write(out+fname,format='MSEED',encoding=11)
    
    
