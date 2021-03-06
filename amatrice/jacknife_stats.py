from numpy import zeros,genfromtxt,mean,std,savetxt,where,isnan
from string import rjust

path='/Users/dmelgar/Slip_inv/Amatrice_M6.6/output/inverse_models/models/'
iter=200
nfaults=1103
mod=zeros((nfaults,iter))
for k in range(iter):
    kread=rjust(str(k),4,'0')
    m=genfromtxt(u'/Users/dmelgar/Slip_inv/Amatrice_M6.6/output/inverse_models/models/jacknife.'+kread+'.inv')
    slip=(m[:,8]**2+m[:,9]**2)**0.5
    mod[:,k]=slip
    
#Write models
mean_fault=m.copy()
std_fault=m.copy()
cv_fault=m.copy()

mean_fault[:,9]=mean(mod,1)
mean_fault[:,8]=0
std_fault[:,9]=std(mod,1)
std_fault[:,8]=0
cv_fault[:,9]=std(mod,1)/mean(mod,1)
cv_fault[:,8]=0
i=where(isnan(cv_fault[:,9])==True)[0]
cv_fault[i,9]=0

fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e\t%10.1f\t%10.1f\t%8.4f\t%8.4e'
savetxt(path+'jacknife.mean.inv',mean_fault,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
savetxt(path+'jacknife.std.inv',std_fault,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
savetxt(path+'jacknife.cv.inv',cv_fault,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')


