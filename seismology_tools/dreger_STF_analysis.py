from scipy.special import gamma,gammaincc
from scipy.integrate import cumtrapz,trapz
from numpy import exp,arange,argmin,zeros,c_,ones
from numpy.linalg import lstsq
from matplotlib import pyplot as plt

tau=arange(0.1,5+1e-10,0.05)
zeta=0.2 ; dt=0.01 ; tmax=100

t=arange(0,tmax+1e-9,dt)
thresh=98.0/100

plt.figure()
trise=zeros(len(tau))
for k in range(len(tau)):
    sdot=(t**zeta)*exp(-t/tau[k])
    s=cumtrapz(sdot,t,initial=0)
    ss=-tau[k]*(t**(zeta))*((t/tau[k])**(-zeta))*gamma(zeta+1)*gammaincc(zeta+1,t/tau[k])
    
    tt=1e-100
    adjust=-tau[k]*(tt**(zeta))*((tt/tau[k])**(-zeta))*gamma(zeta+1)*gammaincc(zeta+1,tt/tau[k])
    ss=ss-adjust
    ss[0]=0
    
    #total slip
    total=ss[-1]
    #When is what percentage reached?
    fraction=ss/total
    i=argmin(abs(fraction-thresh))
    trise[k]=t[i]

    plt.plot(t,sdot)
    plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Slip Rate')

plt.figure()
plt.scatter(tau,trise)
plt.xlabel(r'$\tau$')
plt.ylabel('Rise time')

plt.show()

x=lstsq(c_[tau,ones(len(tau))],trise)
m=x[0][0]
b=x[0][1]

print 't_rise = %.4f*tau+(%.4f)' %(m,b)

#regress for relationship between tau and rise time
