from matplotlib import pyplot as plt
from numpy import genfromtxt,arange

f=genfromtxt('/Users/dmelgar/Slip_inv/Cascadia_gamma/output/forward_models/_wang.grid')
ux=f[:,4]
uy=f[:,3]
uz=f[:,5]
xout=f[:,1]
yout=f[:,2]

plt.figure(figsize=(16,16))
h=(ux**2+uy**2)**0.5
plt.scatter(xout,yout,c=h,lw=0,s=20,vmin=0.0,vmax=0.035)
plt.colorbar()
i=arange(0,len(h),4)
plt.quiver(xout[i],yout[i],ux[i]/h[i],uy[i]/h[i],pivot='mid',linewidths=0.01, edgecolors=('k'),scale=50)
plt.grid()
plt.title('MudPy solution')

plt.figure(figsize=(20,5))

plt.subplot(131)
plt.scatter(xout,yout,c=ux,lw=0,s=20,vmin=-0.05,vmax=0.03)
plt.colorbar()
plt.title('ux')

plt.subplot(132)
plt.scatter(xout,yout,c=uy,lw=0,s=20,vmin=-0.03,vmax=0.03)
plt.colorbar()
plt.title('uy')

plt.subplot(133)
plt.scatter(xout,yout,c=uz,lw=0,s=20,vmin=-0.015,vmax=0.015)
plt.colorbar()
plt.title('uz')

plt.show()