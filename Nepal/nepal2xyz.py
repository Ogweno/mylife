from numpy import genfromtxt,savetxt,zeros,c_

sta=genfromtxt('/Users/dmelgar/Slip_inv/Nepal_fwd/data/station_info/np.sta',usecols=0,dtype='S')
ll=genfromtxt('/Users/dmelgar/Slip_inv/Nepal_fwd/data/station_info/np.sta',usecols=[1,2])
neu=zeros((len(ll),3))
for k in range(len(sta)):
    tmp=genfromtxt('/Users/dmelgar/Slip_inv/Nepal_fwd/output/forward_models/'+sta[k]+'.static.neu')
    neu[k,:]=tmp.transpose()
savetxt('/Users/dmelgar/Slip_inv/Nepal_fwd/output/forward_models/fwdNP.xyz',c_[ll,neu],fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f')
    
