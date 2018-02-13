from numpy import arange,linspace,c_,savetxt
from string import rjust

N=500
x1,x2=-96.2,-93.89
y1,y2=13.21,15.9494
x=linspace(x1,x2,N)
y=linspace(y1,y2,N)


for k in range(len(x)):
    num=rjust(str(k),3,'0')
    print 'rundata.gaugedata.gauges.append([1'+num+','+str(x[k])+','+str(y[k])+',0., 1.e10])'
    
N=5000
x=linspace(x1,x2,N)
y=linspace(y1,y2,N)
savetxt(u'/Users/dmelgar/Chiapas2017/misc/tsun_bathy.txt',c_[x,y],fmt='%.5f')