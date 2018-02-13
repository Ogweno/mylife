from matplotlib import pyplot as plt
from numpy import where,arange,pi,zeros,r_,convolve,array
from mudpy.analysis import subfault_STFs

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

    
    
tr=3.3
ts=1.7
    
ty,sy=yoffe_reg(tr,ts)   
            
                    
k1=134 ; dt1=13
k2=132 ; dt2=18
k3=130 ; dt3=24.5
k4=128 ; dt4=32.5
k5=126 ; dt5=33.5
ymax=1.2
xl=[5,45]
    
epicenter=array([ 84.708,  28.147,  15.   ])
t,M=subfault_STFs('/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/output/inverse_models/models/review_3.3_final2.0001.inv',epicenter,20,15)
plt.figure()    
ax=plt.subplot(511)

ax.plot(t[k1,:],M[k1,:]/(30e9*1e4*1e4),lw=1.5,c='#4682B4')
ax.plot(ty+dt1,sy*(M[k1,:]/(30e9*1e4*1e4)).max(),'k',lw=1.5,c='#8B0000')
ax.set_xlim(xl)
ax.set_yticklabels([])
ax.set_xticklabels([])
smax=(M[k1,:]/(30e9*1e4*1e4)).max()
smax=round(smax*100)/100
ax.set_ylim([0,ymax])
plt.annotate(str(smax), xy=(38,0.75),xytext=(38,0.75),fontsize=14)
plt.annotate('a', xy=(9,0.75),xytext=(9,0.75),fontsize=14)

ax=plt.subplot(512)
ax.plot(t[k2,:],M[k2,:]/(30e9*1e4*1e4),lw=1.5,c='#4682B4')
ax.plot(ty+dt2,sy*(M[k2,:]/(30e9*1e4*1e4)).max(),'k',lw=1.5,c='#8B0000')
ax.set_xlim(xl)
ax.set_yticklabels([])
ax.set_xticklabels([])
#ax.set_xlabel('Time(s)')
smax=(M[k2,:]/(30e9*1e4*1e4)).max()
smax=round(smax*100)/100
ax.set_ylim([0,ymax])
plt.annotate(str(smax), xy=(38,0.75),xytext=(38,0.75),fontsize=14)
plt.annotate('b', xy=(9,0.75),xytext=(9,0.75),fontsize=14)

ax=plt.subplot(513)
ax.plot(t[k3,:],M[k3,:]/(30e9*1e4*1e4),lw=1.5,c='#4682B4')
ax.plot(ty+dt3,sy*(M[k3,:]/(30e9*1e4*1e4)).max(),'k',lw=1.5,c='#8B0000')
ax.set_xlim(xl)
ax.set_ylabel('Slip rate (m/s)')
ax.yaxis.set_label_position("right")
ax.set_yticklabels([])
ax.set_xticklabels([])
smax=(M[k3,:]/(30e9*1e4*1e4)).max()
smax=round(smax*100)/100
ax.set_ylim([0,ymax])
plt.annotate(str(smax), xy=(38,0.75),xytext=(38,0.75),fontsize=14)
plt.annotate('c', xy=(9,0.75),xytext=(9,0.75),fontsize=14)


ax=plt.subplot(514)
ax.plot(t[k4,:],M[k4,:]/(30e9*1e4*1e4),lw=1.5,c='#4682B4')
ax.plot(ty+dt4,sy*(M[k4,:]/(30e9*1e4*1e4)).max(),'k',lw=1.5,c='#8B0000')
ax.set_xlim(xl)
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_xlabel('Time(s)')
smax=(M[k4,:]/(30e9*1e4*1e4)).max()
smax=round(smax*100)/100
ax.set_ylim([0,ymax])
plt.annotate(str(smax), xy=(38,0.75),xytext=(38,0.75),fontsize=14)
plt.annotate('d', xy=(9,0.75),xytext=(9,0.75),fontsize=14)
#ax.set_title('SF '+str(k4))
#ax.legend(['Slip model','Yoffe'],bbox_to_anchor=(2.2, 1),ncol=1,frameon=False,fontsize=14)
    
ax=plt.subplot(515)
ax.plot(t[k5,:]-5,M[k5,:]/(30e9*1e4*1e4),lw=1.5,c='#4682B4')
ax.plot(ty+dt5,sy*(M[k5,:]/(30e9*1e4*1e4)).max(),'k',lw=1.5,c='#8B0000')
ax.set_xlim(xl)
ax.set_yticklabels([])
ax.set_xticklabels(['','10','','20','','30','','40'])
ax.set_xlabel('Time (s) ')
smax=(M[k5,:]/(30e9*1e4*1e4)).max()
smax=round(smax*100)/100
ax.set_ylim([0,ymax])
plt.annotate(str(smax), xy=(38,0.75),xytext=(38,0.75),fontsize=14)
plt.annotate('e', xy=(9,0.75),xytext=(9,0.75),fontsize=14)
ax.legend(['Slip model','Yoffe'],bbox_to_anchor=(1.0, 5.8),ncol=1,frameon=False,fontsize=14)
    
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8, wspace=0.0, hspace=0)
plt.show()