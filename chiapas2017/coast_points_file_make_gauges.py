from numpy import genfromtxt,arange
from scipy.interpolate import interp1d
from string import rjust

dlon=0.01
lon_interp=arange(-97.,-90.,dlon)
delta_lat=0.06

p=genfromtxt('/Users/dmelgar/Chiapas2017/tsunami/coast_points/coast_points')

f=interp1d(p[:,0], p[:,1], kind='linear')
lat_interp=f(lon_interp)-delta_lat

for k in range(len(lon_interp)):
    num=rjust(str(k),3,'0')
    print 'rundata.gaugedata.gauges.append([2'+num+','+str(lon_interp[k])+','+str(lat_interp[k])+',0., 1.e10])'
