from okada_wrapper import dc3dwrapper
from numpy import linspace,zeros,arange,array,cos,sin,deg2rad,c_,savetxt
from matplotlib import pyplot as plt


#Coordinates of observing points in km
x=linspace(-300,300,100)
y=linspace(-300,300,100)

#x=linspace(-50,50,100)
#y=linspace(-50,50,100)

##Fault description
Mw=8.2
M0=10**(Mw*1.5+9.1)
mu=30e9
L=10**(-1.91+0.49*Mw)
W=10**(-1.20+0.36*Mw)
depth=58
area=L*W*1e3*1e3
slip=-M0/(mu*area)
alpha=2./3 #This is alwsays 2/3
strike=310-180
dip=10

#Mw=8.2
#M0=10**(Mw*1.5+9.1)
#mu=30e9
#L=10**(-1.91+0.49*Mw)
#W=10**(-1.20+0.36*Mw)
#depth=58
#area=L*W*1e3*1e3
#ss_slip=0
#ds_slip=-M0/(mu*area)
#alpha=2./3 #This is alwsays 2/3
#strike=310-180
#dip=10

#Mw=7.5
#L=80.0
#W=38.0
#M0=10**(Mw*1.5+9.1)
#mu=30e9
#depth=20.0
#area=L*W*1e3*1e3
#slip=-M0/(mu*area)
#alpha=2./3
#strike=310
#dip=80



#Rotation matrices because code only works  for a north south fault so you need 
#to "rotate" positions by strike angle and then"un rotate" output
theta=strike-90
theta=deg2rad(theta)
R=array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
R2=array([[cos(-theta),-sin(-theta)],[sin(-theta),cos(-theta)]])

#Looping is bad but code is not vectorized
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
                                 [0, slip, 0.0])
        
        # "un-rotate"
        urot=R2.dot(array([[u[0]], [u[1]]]))
        u[0]=urot[0]
        u[1]=urot[1]
        
        #populate output vector
        ux[k]=u[0]
        uy[k]=u[1]
        uz[k]=u[2]
        k+=1
   

        
plt.figure(figsize=(19,7))

plt.subplot(121)
h=(ux**2+uy**2)**0.5
hmax=h.max()
plt.scatter(xout,yout,c=h,lw=0,s=200,vmin=0.0,vmax=hmax,cmap=plt.cm.jet)
plt.colorbar()
i=arange(0,len(h),15)
plt.quiver(xout[i],yout[i],ux[i]/h[i],uy[i]/h[i],pivot='mid',linewidths=0.01, facecolors='w',scale=50)
plt.xlim([x.min(),x.max()])
plt.ylim([y.min(),y.max()])
plt.grid()

plt.subplot(122)
h=uz
hmax=abs(h.max())
plt.scatter(xout,yout,c=h,lw=0,s=200,vmin=-hmax,vmax=hmax,cmap=plt.cm.seismic)
plt.colorbar()
plt.xlim([x.min(),x.max()])
plt.ylim([y.min(),y.max()])
plt.grid()

plt.title('strike=%d' % (strike))
plt.show()

#out=c_[xout,yout,ux,uy,uz]
#savetxt('/Users/dmelgar/GlarmS/Sensitivity/CA_M6.5_2.0cm.txt',out,fmt='%.6f',header='x(km),y(km),ux(m),uy(m),uz(m)')