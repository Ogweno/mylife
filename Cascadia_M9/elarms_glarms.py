from numpy import genfromtxt,arange,savetxt,c_,where
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mudpy import viewFQ

elarms_t=[12.620,13.760,14.240,15.180,16.380,20.520,21.020,21.580,21.860,24.540,28.180,29.020,31.560,32.500,35.440,300]
elarms_m=[7.95,6.97,6.90,6.83,7.14,7.37,7.53,7.53,7.54,7.57,7.58,7.58,7.58,7.58,7.58,7.58]

glarms_t=[30,163]
glarms_m=[8.48,8.72]

##interpolate elarms
t=arange(0,300,0.5)
f = interp1d(elarms_t, elarms_m,bounds_error=False,fill_value=0)
Eout=f(t)
i=where(t<12.66)[0]
Eout[i]=0
savetxt(u'/Users/dmelgar/code/GMT/Cascadia_M9/elarms.mag',c_[t,Eout],fmt='%.2f')
#
#i=where(t<30)[0]
#Eout[i]=0
#i=where(t>=30)[0]
#Eout[i]=8.48
#i=where(t>163)[0]
#Eout[i]=8.72
#savetxt(u'/Users/dmelgar/code/GMT/Cascadia_M9/glarms.mag',c_[t,Eout],fmt='%.2f')
#
#
#glarms=genfromtxt('/Users/dmelgar/Cascadia_M9/glarms_fq.txt')

#t,Mrate=viewFQ.source_time_function(u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.000006.rupt',[-124.616004,45.863800,19.84],dt=3.0)
#ti=arange(0,300,0.5)
#f = interp1d(t, Mrate,bounds_error=False,fill_value=0)
#Mout=f(ti)
#savetxt(u'/Users/dmelgar/code/GMT/Cascadia_M9/stf_inter',c_[ti,(Mout/Mout.max())*2.0+6.0],fmt='%.2f')


#plt.figure()
#plt.plot(elarms_t,elarms_m)
#plt.plot(glarms[:,2],glarms[:,0])
#plt.legend(['ElarmS','GlarmS'])
#plt.xlabel('Seconds since OT')
#plt.ylabel('Magnitude')
#plt.plot([0,300],[8.7,8.7],'--k',lw=2)
#plt.ylim([2,9.5])
#plt.show()