from numpy import genfromtxt,zeros,unique,meshgrid,reshape,arange,c_,savetxt
from scipy.interpolate import interp2d

xl=[20.2,20.85]
yl=[38.3,39.1]
dl=0.001
xnew = arange(xl[0],xl[1], dl)
ynew = arange(yl[0],yl[1], dl)
X,Y=meshgrid(xnew,ynew)


sta=genfromtxt('/Users/dmelgar/Slip_inv/Lefkada_tsun/data/station_info/sfgrid.sta',usecols=0,dtype='S')
lonlat=genfromtxt('/Users/dmelgar/Slip_inv/Lefkada_tsun/data/station_info/sfgrid.sta',usecols=[1,2])
path=u'/Users/dmelgar/Slip_inv/Lefkada_tsun/output/forward_models/'
fout=u'/Users/dmelgar/Slip_inv/Lefkada_tsun/output/forward_models/_3d_deformation.txt'

n=zeros(len(sta))
e=zeros(len(sta))
z=zeros(len(sta))

for k in range(len(sta)):
    neu=genfromtxt(path+sta[k]+'.static.neu')
    n[k]=neu[0]
    e[k]=neu[1]
    z[k]=neu[2]
    
#interpolation
x=unique(lonlat[:,0])
y=unique(lonlat[:,1])
xx, yy = meshgrid(x, y)
zn = reshape(n,(len(x),len(y))).T
ze = reshape(e,(len(x),len(y))).T
zz = reshape(z,(len(x),len(y))).T
fn = interp2d(x, y, zn, kind='linear',bounds_error=False,fill_value=0)
fe = interp2d(x, y, ze, kind='linear',bounds_error=False,fill_value=0)
fz = interp2d(x, y, zz, kind='linear',bounds_error=False,fill_value=0)

nout=fn(xnew,ynew)
eout=fe(xnew,ynew)
zout=fz(xnew,ynew)

out=c_[X.ravel(),Y.ravel(),nout.ravel(),eout.ravel(),zout.ravel()]
savetxt(fout,out,fmt='%.4f',header='lon,lat,N(m),E(m),Up(m)')