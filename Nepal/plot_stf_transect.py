from mudpy.analysis import subfault_STFs
from numpy import array,arange,zeros,where,pi,r_,arange
from scipy.signal import convolve
from matplotlib import pyplot as plt

#rupture=u'/Users/dmelgar/Slip_inv/Nepal_ttests_inv/output/inverse_models/models/4s_vr3.2.0000.inv'
rupture='/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/output/inverse_models/models/review_3.2.0000.inv'
epicenter=array([ 84.708,  28.147,  15.   ])

t,M=subfault_STFs(rupture,epicenter,20,15)

def H(t):
    unit_step = arange(t.shape[0])
    lcv = arange(t.shape[0])
    for place in lcv:
        if t[place] == 0:
           unit_step[place] = .5
        elif t[place] > 0:
           unit_step[place] = 1
        elif t[place] < 0:
            unit_step[place] = 0
    return unit_step

def yoffe_reg(tr,ts):
    dt=0.01
    t=arange(dt,100,dt)
    #Make Joffe function
    H1=H(t)
    H2=H(tr-t)
    Y=(2./(pi*tr))*H1*H2*abs(((tr-t)/t))**0.5

    #Triangle
    W=(1./ts**2)*(t*H(t)*H(ts-t)+(2*ts-t)*H(t-ts)*H(2*ts-t))
    i=where(t==ts)[0]
    W[i]=max(W)

    #Convolve
    s=convolve(Y,W)
    s=(s[0:len(t)]/s.max())
    s=r_[zeros(60/dt),s]
    t=arange(-60,len(s)*dt-60,dt)
    return t,s

vr=3.2
trise=4
tr=3.3
ts=1.7
xl=[0,60]
ty,sy=yoffe_reg(tr,ts)  
ymax=1.3

subfaults=range(123,139,2)
#dt=arange(0,len(subfaults))*(18./vr)+3.14
#dt=[3.14,7.02,12.99,19.1,25.1,31.2,37.3,43.8]  #2.8
#dt=[3.14,7.02,12.99,19.1,25.1,31.2,37.3,43.8]  #3.0
dt=[3.14,7.02,12.99,19.1,25.1,31.2,37.3,43.8]  #3.2
#dt=[3.14,7.02,12.99,19.1,25.1,31.2,37.3,43.8]  #3.4
#dt=[3.14,7.02,12.99,19.1,25.1,31.2,37.3,43.8]  #3.6
dt=dt[::-1]
fig, axarr = plt.subplots(len(subfaults), 1) 
k=0
for sub in subfaults:
    ax=axarr[k]
    smax=(M[sub,:]/(30e9*1e4*1e4)).max()
    smax=round(smax*100)/100
    ax.plot(t[sub,:],M[sub,:]/(30e9*1e4*1e4)/smax,lw=1.5,c='#4682B4')
    ax.plot(ty+dt[k],sy*(M[sub,:]/(30e9*1e4*1e4)).max()/smax,'k',lw=1.5,c='#8B0000')
    ax.set_xlim(xl)
    ax.set_yticklabels([])
    if k<(len(subfaults)-1):
        ax.xaxis.set_ticklabels([])
    if k==(len(subfaults)-1):
        ax.set_xlabel('Seconds after OT')

    ax.set_ylim([0,1.1])
    ax.annotate('sub'+str(sub), xy=(50,0.2),fontsize=12)  
    ax.annotate(str(smax), xy=(50,0.8),fontsize=12)  
    k+=1

axarr[0].set_title(r'STFs (m/s), $\tau_r=$'+str(trise)+', $v_r=$'+str(vr))
plt.show()
    