from numpy import genfromtxt,savetxt,arange,r_

f=genfromtxt('/Users/dmelgar/Slip_inv/Amatrice_3plane/data/model_info/south_segment.fault')
i=arange(0,11)
out='/Users/dmelgar/Slip_inv/Amatrice_3plane/data/model_info/south_segment_prune.fault'

for k in range(12):
    if k==0:
        fout=f[i,:]
        i=i+20
    else:
        fout=r_[fout,f[i,:]]
        i=i+20
        
savetxt(out,fout,fmt='%d\t%.6f\t%.6f\t%.3f\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f')