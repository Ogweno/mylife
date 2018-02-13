from numpy import linspace,savetxt,c_

x1,y1 = -95.2738,13.4614
x2,y2 = -93.60,15.9789

npts=500

x1=360+x1
x2=360+x2
xout=linspace(x1,x2,npts)
yout=linspace(y1,y2,npts)

savetxt('/Users/dmelgar/Chiapas2017/misc/AA360.xy',c_[xout,yout],fmt='%.6f')