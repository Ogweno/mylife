from mudpy.analysis import MT
from numpy import genfromtxt,zeros,where
from obspy.imaging.scripts import mopad
from matplotlib import pyplot as plt

mts=genfromtxt(u'/Users/dmelgar/Chiapas2017/seismicity/outer_rise.psmeca')
st=zeros(len(mts)*2)
di=zeros(len(mts)*2)
ra=zeros(len(mts)*2)
for k in range(len(mts)):
    depth=mts[k,2]
    mt=mts[k,3:9]
    MT=mopad.MomentTensor(mt,system='USE')
    np1,np2=MT.get_fps()
    #mt=MT(mrr,mrt,mrf,mtt,mtf,mff,lon,lat,depth,mt_style='xyz')
    #mt.get_nodal_planes()
    st[2*k]=np1[0]
    st[2*k+1]=np2[0]
    di[2*k]=np1[1]
    di[2*k+1]=np2[1]
    ra[2*k]=np1[2]
    ra[2*k+1]=np2[2]
    
i=where((ra>-135)&(ra<-45))[0]



mts=genfromtxt(u'/Users/dmelgar/Chiapas2017/seismicity/not_outer_rise.psmeca')
st2=zeros(len(mts)*2)
di2=zeros(len(mts)*2)
ra2=zeros(len(mts)*2)
for k in range(len(mts)):
    depth=mts[k,2]
    mt=mts[k,3:9]
    MT=mopad.MomentTensor(mt,system='USE')
    np1,np2=MT.get_fps()
    #mt=MT(mrr,mrt,mrf,mtt,mtf,mff,lon,lat,depth,mt_style='xyz')
    #mt.get_nodal_planes()
    st2[2*k]=np1[0]
    st2[2*k+1]=np2[0]
    di2[2*k]=np1[1]
    di2[2*k+1]=np2[1]
    ra2[2*k]=np1[2]
    ra2[2*k+1]=np2[2]
    
i2=where((ra2>-135)&(ra2<-45))[0]


plt.figure(figsize=(16,8))

plt.subplot(221)
plt.hist(st[i],15,color='#3399FF')
plt.xlim([0,360])
plt.legend(['Outer rise'],loc=2)
plt.ylabel('Frecuency')
plt.subplot(222)
plt.hist(di[i],15,color='#3399FF')
plt.xlim([0,90])




plt.subplot(223)
plt.hist(st2[i2],15,color='#FFB266')
plt.xlim([0,360])
plt.xlabel('Strike')
plt.ylabel('Frecuency')
plt.legend(['In slab'],loc=2)
plt.plot([309,309],[0,96],'--',lw=2)

plt.subplot(224)
plt.hist(di2[i2],15,color='#FFB266')
plt.xlabel('Dip')
plt.xlim([0,90])
plt.plot([77,77],[0,38],'--',lw=2)

plt.show()