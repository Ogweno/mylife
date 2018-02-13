from numpy import linspace,savetxt,c_

x1,y1 = -96.338,14.036
x2,y2 = -95.84,16.51

npts=500

x1=x1
x2=x2
xout=linspace(x1,x2,npts)
yout=linspace(y1,y2,npts)

savetxt('/Users/dmelgar/Chiapas2017/misc/BB.xy',c_[xout,yout],fmt='%.6f')