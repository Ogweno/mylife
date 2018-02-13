from numpy import genfromtxt,zeros,sqrt,argmin,c_,savetxt,arange,r_
from mudpy import forward

strike=163
dip=74
rake=180

f=genfromtxt('/Users/dmelgar/Napa2014/stochastic/inversion/v_fine1.0014.inv.total')
fout='/Users/dmelgar/Napa2014/stochastic/inversion/napa_200m.fault'
fout_rupt='/Users/dmelgar/Napa2014/stochastic/inversion/napa_200m.rupt'



lonlat=f[:,1:3]
ss=f[:,8]
ds=f[:,9]
S=(ss**2+ds**2)**0.5

#now make geometry of finely discretized fault
forward.makefault(fout,strike,dip,173,0.2,0.2,[-122.3174,38.2118,10.729],55,22,0.5)
f=genfromtxt(fout)

#Make mean rupture model
slip_out=zeros((len(f),1))
for k in range(len(f)):
    if k%100==0:
        print k
    dist=sqrt((f[k,1]-lonlat[:,0])**2+(f[k,2]-lonlat[:,1])**2)
    slip_out[k,0]=S[argmin(dist)]
    

#trim edges
i=arange(50,173)
iout=i.copy()
for k in range(1,78):
    iout=r_[iout,i+(173*k)]


out=c_[arange(len(iout))+1,f[iout,1:8],slip_out[iout],zeros((len(iout),1)),f[iout,-2:],zeros((len(iout),2))]
#out=c_[f[:,0:8],slip_out,zeros((len(f),1)),f[:,-2:],zeros((len(f),2))]
savetxt(fout_rupt,out,fmt='%10d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t')
