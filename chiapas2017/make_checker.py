from numpy import genfromtxt,savetxt,unique,array,r_

fout='/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/forward_models/checker_new.rupt'
f=genfromtxt(u'/Users/dmelgar/Slip_inv/Chiapas_hernandez_new/output/inverse_models/models/_final/final_vr_3.4.0007.inv')

#Keep unique subfaults
i=unique(f[:,0])
f=f[0:len(i),:]
f[:,8]=0
f[:,9]=0

#Arrange
i1=array([0,1,4,5,8,9,12,13,16,17])
i2=i1+19
i=r_[i1,i2]
f[i,9]=-1

i1=i2+21
i1=i1[:-1]
i2=i1+19
i=r_[i1,i2]
f[i,9]=1

i1=i2+17
i1=r_[i1,i1[-1]+1]
i2=i1+19
i=r_[i1,i2]
f[i,9]=1

i1=i2+21
i1=i1[0:-2]
i1=r_[i1,i1[-1]+3]
i2=i1+19
i=r_[i1,i2]
f[i,9]=1

fmt='%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%.2f\t%.2f\t%.4f\t%.4e'
savetxt(fout,f,fmt=fmt)