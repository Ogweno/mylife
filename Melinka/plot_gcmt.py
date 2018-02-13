from numpy import genfromtxt,log10
from matplotlib import pyplot as plt

cat=genfromtxt('/Users/dmelgar/code/GMT/Melinka/gcmt_thrust.psmeca')

# mrr, mtt, mff, mrt, mrf, mtf, exp
exp=10**cat[:,9]

mrr=cat[:,3]*exp
mtt=cat[:,4]*exp
mff=cat[:,5]*exp
mrt=cat[:,6]*exp
mrf=cat[:,7]*exp
mtf=cat[:,8]*exp

mtr=mrt
mfr=mrf
mft=mtf


M0=(1./(2**0.5))*((mrr**2+mtt**2+mff**2+mrt**2+mrf**2+mtf**2+mtr**2+mfr**2+mft**2)**0.5)/1e7
Mw=(2./3)*(log10(M0)-9.1)

lon=cat[:,0]
lat=cat[:,1]

plt.stem(lat,Mw,'k')
plt.stem([-43.416,-43.416],[7.6,7.6],'r')
plt.ylim([5,8])
plt.xlim([-44,-36])
plt.ylabel(r'$M_w$')
plt.xlabel('Latitude')

plt.show()