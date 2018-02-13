from numpy import genfromtxt,savetxt,unique,array,r_

fout='/Users/dmelgar/Slip_inv/Melinka_usgs/forward_models/checker.rupt'
f=genfromtxt(u'/Users/dmelgar/Slip_inv/Melinka_usgs/output/inverse_models/models/gps_sm_insar_w3.01.0001.inv')

#Keep unique sunfaults
i=unique(f[:,0])
f=f[0:len(i),:]
f[:,8]=0
f[:,9]=0

#Arrange
i1=array([0,1,2,6,7,8,12,13,14,18,19,20])
i2=i1+24
i3=i2+24
i=r_[i1,i2,i3]
f[i,9]=1

i1=i3+27
i2=i1+24
i3=i2+24
i=r_[i1,i2,i3]
f[i,9]=1

i1=i3+21
i2=i1+24
i3=i2+24
i=r_[i1,i2,i3]
f[i,9]=1

i1=i3+27
i2=i1+24
i3=i2+24
i=r_[i1,i2,i3]
f[i,9]=1

i1=i3+21
i2=i1+24
i3=i2+24
i=r_[i1,i2,i3]
f[i,9]=1

fmt='%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%.2f\t%.2f\t%.4f\t%.4e'
savetxt(fout,f,fmt=fmt)