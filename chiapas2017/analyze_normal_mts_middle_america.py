from mudpy.analysis import MT
from numpy import genfromtxt,zeros,where
from obspy.imaging.scripts import mopad
from matplotlib import pyplot as plt

mts=genfromtxt(u'/Users/dmelgar/Chiapas2017/seismicity/normal_faults_middle_america.txt')
st=zeros(len(mts)*2)
di=zeros(len(mts)*2)
ra=zeros(len(mts)*2)
z=zeros(len(mts)*2)
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
    z[k]=depth
    
i=where((ra>-135)&(ra<-45))[0]
di=di[i]
z=z[i]

bins=10
plt.figure()

i=where((z<20)&(z>0))[0]
plt.hist(di[i])

i=where((z<150)&(z>20))[0]
plt.hist(di[i],histtype='step')

i=where((z<150)&(z>20))[0]
plt.hist(di[i],histtype='step')

plt.legend(['shallow','inter','deep'])
plt.show()