from numpy import zeros,genfromtxt,mean,std,savetxt,where,isnan,unique
from string import rjust
from mudpy import gmttools

path='/Users/dmelgar/Slip_inv/puebla_herloc_NP2/output/inverse_models/models/'
iter=86
nfaults=27*17
mod=zeros((nfaults,iter))
for k in range(iter):
    kread=rjust(str(k),4,'0')
    gmttools.make_total_model(u'/Users/dmelgar/Slip_inv/puebla_herloc_NP2/output/inverse_models/models/jacknife.'+kread+'.inv',0)
    m=genfromtxt(u'/Users/dmelgar/Slip_inv/puebla_herloc_NP2/output/inverse_models/models/jacknife.'+kread+'.inv.total')
    slip=(m[:,8]**2+m[:,9]**2)**0.5
    mod[:,k]=slip
    
#Write models
mean_fault=m[0:nfaults,:].copy()
std_fault=m[0:nfaults,:].copy()
cv_fault=m[0:nfaults,:].copy()

mean_fault[:,9]=mean(mod,1)
mean_fault[:,8]=0
std_fault[:,9]=std(mod,1)
std_fault[:,8]=0
cv_fault[:,9]=std(mod,1)/mean(mod,1)
cv_fault[:,8]=0
i=where(isnan(cv_fault[:,9])==True)[0]
cv_fault[i,9]=0

fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e\t%10.1f\t%10.1f\t%8.4e'
savetxt(path+'jacknife.mean.inv.total',mean_fault,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
savetxt(path+'jacknife.std.inv.total',std_fault,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
savetxt(path+'jacknife.cv.inv.total',cv_fault,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')


