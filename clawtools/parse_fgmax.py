import sys
from numpy import genfromtxt,unique,where,zeros,c_,savetxt,nan,meshgrid,arange,double,isnan,log10
from scipy.interpolate import griddata
from matplotlib import pyplot as plt
from matplotlib import cm


#fgmax_file='/Users/dmelgar/Tsunamis/Cascadia/fort.FG2.valuemax'
#aux_file='/Users/dmelgar/Tsunamis/Cascadia/fort.FG2.aux1'
#out_file='/Users/dmelgar/Tsunamis/Cascadia/etamax'

#fgmax_file=u'/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/fort.FG1.valuemax'
#aux_file=u'/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/fort.FG1.aux1'
#out_file='/Users/dmelgar/NewZealand2016/tsunami/max_eta1.txt'

fgmax_file=u'/Users/dmelgar/Tsunamis/tehuantepec/_output/fort.FG1.valuemax'
aux_file=u'/Users/dmelgar/Tsunamis/tehuantepec/_output/fort.FG1.aux1'
out_file='/Users/dmelgar/Tsunamis/tehuantepec/max_eta.txt'

outlog=False
amr=1

#fgmax_file=sys.argv[1]
#aux_file=sys.argv[2]
#out_file=sys.argv[3]

#Anything larger than this is clipped
maxeta=10

print 'Reading fgmax data from '+fgmax_file
print 'Reading auxiliary data from '+aux_file
print 'Will save output to '+out_file
print 'Clipping eta at '+str(maxeta)+'m'

fgmax=genfromtxt(fgmax_file)
aux=genfromtxt(aux_file)

#Get coordinates
lon=fgmax[:,0]
lat=fgmax[:,1]

#Maximum elevation
H=fgmax[:,3]

#At which grid level?
amr_level=fgmax[:,2]

#Make eta from info at appropriate grid level
num_amr=unique(amr_level)
eta=zeros(len(H))
b=zeros(len(H))
for k in range(len(num_amr)):
    i=where(amr_level==num_amr[k])[0]
    b[i]=aux[i,int(1+num_amr[k])]# Use data at that amr level
#Add em up
eta=H+b 

#Filter by AMR
i=where(amr_level==amr)[0]
lon=lon[i]
lat=lat[i]
eta=eta[i]
b=b[i]
#    
#Set dry cells to NaN and clip at etamax level
i=where(b>0)[0]
eta[i]=nan
#i=where(eta>maxeta)[0]
#eta[i]=nan


#Grid tor egular
delta=0.005
xi=arange(lon.min(),lon.max()+0.00001,delta)
yi=arange(lat.min(),lat.max()+0.00001,delta)
X,Y=meshgrid(xi,yi)
Z=griddata((lon,lat),eta,(X,Y))
#Write to file
xout=zeros(Z.size)
yout=zeros(Z.size)
zout=zeros(Z.size)
k=0
for kx in range(xi.size):
    for ky in range(yi.size):
        xout[k]=X[ky,kx]
        yout[k]=Y[ky,kx]
        zout[k]=Z[ky,kx]
        k+=1
#i=where(isnan(zout)==1)[0]
i=arange(len(zout))
#zout[i]=0
if outlog:
    out=c_[xout[i],yout[i],log10(zout[i])]
else:
    out=c_[xout[i],yout[i],zout[i]]
#if outlog:
#    out=c_[xout,yout,log10(zout)]
#else:
#    out=c_[xout,yout,zout]

    
    
    
savetxt(out_file,out,fmt='%.8f\t%.8f\t%.6f',header='lon,lat,eta(m)')


#Plot to check you're not dumb
plt.figure()
plt.scatter(xout,yout,c=zout,cmap=cm.rainbow,marker='s',s=50,vmin=0,vmax=2,linewidth=0)
plt.colorbar()
plt.show()