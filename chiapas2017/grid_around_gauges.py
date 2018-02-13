from numpy import arange,meshgrid
from string import rjust

f=open('/Users/dmelgar/DEMs/SRTM15/grid_search.gauge','w')


#puerto angel
dl=1./3600
N=12
center=[-96.493,15.6647]

x=arange(center[0]-N/2*dl,center[0]+N/2*dl,dl)
y=arange(center[1]-N*dl,center[1]+N/2*dl,dl)
X,Y=meshgrid(x,y)
x=X.ravel()
y=Y.ravel()

for k in range(len(x)):
    num=rjust(str(k),4,'0')
    f.write('rundata.gaugedata.gauges.append([1'+num+','+str(x[k])+','+str(y[k])+',0., 1.e10])\n')




#hautulco
dl=1./3600
N=12
center=[-96.130000,15.753333]

x=arange(center[0],center[0]+N*dl,dl)
y=arange(center[1]-N*dl,center[1],dl)
X,Y=meshgrid(x,y)
x=X.ravel()
y=Y.ravel()

for k in range(len(x)):
    num=rjust(str(k),4,'0')
    f.write('rundata.gaugedata.gauges.append([2'+num+','+str(x[k])+','+str(y[k])+',0., 1.e10])\n')





#salina
dl=1./3600
N=12
center=[-95.1964,16.1676]

x=arange(center[0]-N/2*dl,center[0]+N*dl,dl)
y=arange(center[1]-N*dl,center[1]+N/2*dl,dl)
X,Y=meshgrid(x,y)
x=X.ravel()
y=Y.ravel()

for k in range(len(x)):
    num=rjust(str(k),4,'0')
    f.write('rundata.gaugedata.gauges.append([3'+num+','+str(x[k])+','+str(y[k])+',0., 1.e10])\n')





#puerto madero
dl=1./3600
N=12
center=[-92.403,14.711]

x=arange(center[0]-N/2*dl,center[0]+N/2*dl,dl)
y=arange(center[1]-N/2*dl,center[1]+N/2*dl,dl)
X,Y=meshgrid(x,y)
x=X.ravel()
y=Y.ravel()

for k in range(len(x)):
    num=rjust(str(k),4,'0')
    f.write('rundata.gaugedata.gauges.append([4'+num+','+str(x[k])+','+str(y[k])+',0., 1.e10])\n')




# puerto madero 2
dl=1./3600
N=12
center=[-92.4081,14.7004]

x=arange(center[0]-N/2*dl,center[0]+N/2*dl,dl)
y=arange(center[1]-N/2*dl,center[1]+N/2*dl,dl)
X,Y=meshgrid(x,y)
x=X.ravel()
y=Y.ravel()

for k in range(len(x)):
    num=rjust(str(k),4,'0')
    f.write('rundata.gaugedata.gauges.append([5'+num+','+str(x[k])+','+str(y[k])+',0., 1.e10])\n')



f.close()