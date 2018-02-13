from numpy import genfromtxt,where,linspace
from scipy.interpolate import interp1d

f=genfromtxt('/Users/dmelgar/DART_analysis/first_wave.txt')

i=where((f[:,0]>-110) & (f[:,0]<-86) & (f[:,1]>10) & (f[:,1]<24))[0]

f=f[i,:]
x=linspace(-110,-86,300)
I=interp1d(f[:,1],f[:,2]