from numpy import genfromtxt,ones,where,savetxt,arange
from pyproj import Geod
from mudpy import forward

dlim=490.0 #in km

#Read fault model
fault=genfromtxt(u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.000006.rupt')
fout=u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.short.000006.rupt'
log=u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.000006.log'

#Look for rise time=0 faults and get rid of them then re-write consecutive subfault numbers
i=where(fault[:,7]>0)[0]
fault=fault[i,:]
fault[:,0]=arange(len(fault))

savetxt(fout,fault,fmt='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e',header='#')


#Make srf
forward.mudpy2srf(fout,log)
