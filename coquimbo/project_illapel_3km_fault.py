from numpy import genfromtxt,sqrt,zeros,argmin,c_,savetxt,where,arange
from scipy.interpolate import Rbf
from matplotlib import pyplot as plt

rupt=genfromtxt('/Users/dmelgar/Slip_inv/Coquimbo_4s/output/inverse_models/models/gps_sm_tg_insar_9win_vel2.0_final2.0008.inv.total')
f=genfromtxt('/Users/dmelgar/Coquimbo2015/coquimbo_FQ.fault')
fout_fault='/Users/dmelgar/Coquimbo2015/coquimbo_FQ_3km.fault'
fout_rupt='/Users/dmelgar/Coquimbo2015/coquimbo_FQ_3km.rupt'


#P=Rbf(rupt[:,2],rupt[:,1],S)
#slip=P(f[:,1],f[:,2])
#
slip=zeros((len(f),1))
S=sqrt(rupt[:,8]**2+rupt[:,9]**2)
for k in range(len(f)):
    if k%100==0:
        print k
    dist=sqrt((f[k,1]-rupt[:,1])**2+(f[k,2]-rupt[:,2])**2)
    slip[k,0]=S[argmin(dist)]
    
#min slip filter
i=where(slip>0.5)[0]
f=f[i,:]
f[:,0]=arange(1,len(i)+1,1)
slip=slip[i]

plt.scatter(f[:,1],f[:,2],c=slip,lw=0,s=5,cmap=plt.cm.nipy_spectral_r)
plt.show()

savetxt(fout_fault,f,fmt='%10d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t')

out=c_[f[:,0:8],zeros((len(f),1)),slip,f[:,-2:],zeros((len(f),2))]
savetxt(fout_rupt,out,fmt='%10d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t')

#savetxt('/Users/dmelgar/Tohoku2011/Minsons/minson_6km.fault',f,fmt='%10d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f')