from numpy import genfromtxt,zeros,array,log10,where,savetxt
from numpy.linalg import norm

mt=genfromtxt('/Users/dmelgar/Coquimbo2015/afters/MTs.txt')
Mw=zeros(len(mt))
for k in range(len(mt)):
    m=array([[mt[k,3],mt[k,6],mt[k,7]],[mt[k,6],mt[k,4],mt[k,8]],[mt[k,7],mt[k,8],mt[k,8]]])*(10**mt[k,9])
    M0=norm(m,'fro')/2**0.5
    Mw[k]=(2./3)*(log10(M0)-9.1)
i=where(Mw>=6)[0]
mtout=mt[i,:]
savetxt('/Users/dmelgar/Coquimbo2015/afters/large_MTs.txt',mtout,fmt='%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%4i')