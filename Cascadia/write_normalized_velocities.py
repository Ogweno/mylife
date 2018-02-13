from numpy import genfromtxt,arange,savetxt,c_

#f=genfromtxt('/Users/dmelgar/Slip_inv/Cascadia_gauss/output/forward_models/_gauss.grid')
#fout=u'/Users/dmelgar/code/GMT/Cascadia_seafloor/gauss_vectors.txt'
f=genfromtxt('/Users/dmelgar/Slip_inv/Cascadia_gamma/output/forward_models/_wang.grid')
fout=u'/Users/dmelgar/code/GMT/Cascadia_seafloor/wang_vectors.txt'
factor=11



#normalize and decimate
i=arange(0,len(f),factor)
ux=f[i,4]
uy=f[i,3]
uz=f[i,5]
xout=f[i,1]
yout=f[i,2]
h=(ux**2+uy**2)**0.5
ux=ux/h
uy=uy/h

#write file
out=c_[xout,yout,ux,uy]
savetxt(fout,out,fmt='%.6f\t%.6f\t%.6f\t%.6f')
