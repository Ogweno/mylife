from numpy import genfromtxt,zeros,sqrt,argmin,c_,savetxt,fliplr
from mudpy import forward

strike=357
dip=18
rake=90

f=genfromtxt('/Users/dmelgar/Iquique2014/GFZ_model/kinematic_0401.dat')
fout='/Users/dmelgar/Iquique2014/GFZ_model/iquique_gfz.fault'
fout_rupt='/Users/dmelgar/Iquique2014/GFZ_model/iquique_gfz.rupt'

n=len(f)/2
#split in two
f1=f[0:n,:]
f2=f[n:,:]

lonlat=f1[:,0:2]
lonlat=fliplr(lonlat)
ss=f1[:,8]
ds=f2[:,8]
S=(ss**2+ds**2)**0.5

#now make geometry of finely discretized fault
#forward.makefault(fout,strike,dip,45,3.0,3.0,[-70.593,-19.70,29.36],13,13,0.5)
forward.makefault(fout,strike,dip,65,3.0,3.0,[-70.593,-19.70,29.36],30,30,0.5)
f=genfromtxt(fout)

#Make mean rupture model
slip_out=zeros((len(f),1))
for k in range(len(f)):
    if k%100==0:
        print k
    dist=sqrt((f[k,1]-lonlat[:,0])**2+(f[k,2]-lonlat[:,1])**2)
    slip_out[k,0]=S[argmin(dist)]
    

out=c_[f[:,0:8],zeros((len(f),1)),slip_out,f[:,-2:],zeros((len(f),2))]
savetxt(fout_rupt,out,fmt='%10d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t')
