from okada_wrapper import dc3dwrapper
from numpy import linspace,zeros,arange,array,cos,sin,deg2rad,r_,expand_dims,c_,savetxt
from matplotlib import pyplot as plt

#fout='/Users/dmelgar/GlarmS/Sensitivity/Cascadia_M6.5_2.0cm.txt'
#fout='/Users/dmelgar/GlarmS/Sensitivity/Cascadia_M7.0_2.0cm.txt'
#fout='/Users/dmelgar/GlarmS/Sensitivity/Cascadia_M7.5_2.0cm.txt'

x=linspace(-10,10,500)
y=linspace(-10,10,500)

#Mw=6.5
#L=21.0
#W=13.0

#Mw=7.0
#L=41.0
#W=22.0

Mw=4
L=0.001
W=0.001

M0=10**(Mw*1.5+9.1)
mu=30e9
depth=6.0
area=L*W*1e3*1e3
slip=M0/(mu*area)*-1
alpha=2./3
strike=330
dip=45




theta=strike-90
theta=deg2rad(theta)
R=array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
R2=array([[cos(-theta),-sin(-theta)],[sin(-theta),cos(-theta)]])

xout=zeros(len(x)*len(y)) 
yout=zeros(len(x)*len(y))
ux=zeros(len(x)*len(y))
uy=zeros(len(x)*len(y))
uz=zeros(len(x)*len(y)) 
k=0                            
for kx in range(len(x)):
    print kx
    for ky in range(len(y)):
        
        xout[k]=x[kx]
        yout[k]=y[ky]
        
        #Calculate on rotated position
        xy=R.dot(array([[x[kx]], [y[ky]]]))
        success, u, grad_u = dc3dwrapper(alpha, [xy[0], xy[1], 0.0],
                                 depth, dip, [-L/2, L/2], [-W/2, W/2],
                                 [0.0, slip, 0.0])
        
        
        urot=R2.dot(array([[u[0]], [u[1]]]))
        u[0]=urot[0]
        u[1]=urot[1]
        
        ux[k]=u[0]
        uy[k]=u[1]
        uz[k]=u[2]
        k+=1
   

        
#plt.figure(figsize=(16,16))
#h=(ux**2+uy**2)**0.5
#plt.scatter(xout,yout,c=h,lw=0,s=200,vmin=0.0,vmax=hmax)
#plt.colorbar()
#i=arange(0,len(h),7)
#plt.quiver(xout[i],youti],ux[i]/h[i],uy[i]/h[i],pivot='mid',linewidths=0.01, edgecolors=('k'),scale=50)
#plt.xlim([x.min(),x.max()])
#plt.ylim([y.min(),y.max()])
#plt.grid()
#plt.title('strike=%d' % (strike))


plt.figure(figsize=(20,5))
plt.subplot(131)
hmax=abs(ux).max()
plt.scatter(xout,yout,c=ux,lw=0,s=20,vmin=-hmax,vmax=hmax)
plt.colorbar()
plt.xlim([x.min(),x.max()])
plt.ylim([y.min(),y.max()])
plt.grid()
plt.title('ux')

plt.subplot(132)
hmax=abs(uy).max()
plt.scatter(xout,yout,c=uy,lw=0,s=20,vmin=-hmax,vmax=hmax)
plt.colorbar()
plt.xlim([x.min(),x.max()])
plt.ylim([y.min(),y.max()])
plt.grid()
plt.title('uy')

plt.subplot(133)
vmax=abs(uz).max()
plt.scatter(xout,yout,c=uz,lw=0,s=20,vmin=-vmax,vmax=vmax)
plt.colorbar()
plt.xlim([x.min(),x.max()])
plt.ylim([y.min(),y.max()])
plt.grid()
plt.title('uz')


plt.show()

