from numpy import genfromtxt,arange,r_,savetxt

f=genfromtxt('/Users/dmelgar/Chiapas2017/slab/chiapas_normal_ssn.fault')

i=arange(11,40)
iout=i.copy()

for k in range(7):
    iout=r_[iout,i+40*(k+1)]
    
fout=f[iout,:]
fout[:,0]=arange(1,len(fout)+1)

fmt='%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f'
savetxt('/Users/dmelgar/Chiapas2017/slab/chiapas_normal_ssn_trim.fault',fout,fmt=fmt)