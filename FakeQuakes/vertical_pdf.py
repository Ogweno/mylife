from glob import glob
from numpy import genfromtxt,zeros,diff,linspace,histogram,where,log
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

folders=glob('/Users/dmelgar/FakeQuakes/Cascadia_final1/output/dreger_4tau/cascadia*')

for k in range(len(folders)):
    f=glob(folders[k]+'/*.offsets')
    print f[0]
    up=genfromtxt(f[0],usecols=5)
    if k==0:
        lonlat=genfromtxt(f[0],usecols=[1,2])
        U=zeros((len(up),len(folders)))
    U[:,k]=up
    
plt.figure()


for k in range(len(U)):
    xinterp=linspace(-1,1,200)
    n,bins=histogram(U[k,:],bins=linspace(-2,2,100))
    bin_center=bins[0:-1]+diff(bins)
    I=interp1d(bin_center,n,fill_value="extrapolate")
    yinterp=I(xinterp)
    yinterp/=yinterp.max()
    yinterp=log(yinterp)
    yinterp*=0.5
    
    
    #Move by lat lon
    xplot=xinterp+lonlat[k,0]
    yplot=yinterp+lonlat[k,1]
    
    #Plot
    i=where(xinterp<0)[0]
    plt.plot(xplot[i],yplot[i],'b',lw=2)
    i=where(xinterp>=0)[0]
    plt.plot(xplot[i],yplot[i],'r',lw=2)
    