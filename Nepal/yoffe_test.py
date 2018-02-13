from numpy import arange,pi,where,convolve,r_,diff,zeros
from matplotlib import pyplot as plt
from obspy import read
import matplotlib
from scipy.integrate import cumtrapz

#read in impulse response
st=read(u'/Users/dmelgar/Slip_inv/Nepal_Yoffe/GFs/dynamic/avouac.mod_13.0920.sub0001/KKN4.subfault0001.DS.vel.zi')
st[0].data=st[0].data/100
dt=0.01
fc=2.5
tr_s=1.0
ts_s=0.2
tr_l=3.3
ts_l=1.7

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

def low_pass_filter(tr,fcorner,dt,order):
    from scipy.signal import filtfilt,butter
    fnyquist=1./(2*dt)
    print fnyquist
    Fc=fcorner/fnyquist
    print Fc
    b, a = butter(order, Fc,btype='lowpass')
    y = filtfilt(b, a, tr)
    return y

#Do short Yoffe 
t,s=yoffe_reg(tr_s,ts_s)
i=where(t>=0)[0]
t=t[i]
s=s[i]
s[0]=0
sd=cumtrapz(s,t,initial=0)
#Convolve
kkn4=convolve(s,st[0].data)*0.031
tout=arange(0,len(kkn4)*dt,dt)
kkn4_d=cumtrapz(kkn4,tout,initial=0)
kkn4=low_pass_filter(kkn4,fc,dt,2)
kkn4_d=low_pass_filter(kkn4_d,fc,dt,2)



fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.plot(tout-0.3, kkn4,'b-')
ax1.set_xlim([-0.3,4])
ax1.set_ylabel('Vel. or slip rate (m/s)')
ax1.plot(t, s, 'k')
ax1.legend(['KKN4','Yoffe'])
ax1.grid()
plt.annotate(r'$\tau_r = '+str(tr_s)+'s$', xy=(2,0.4),xytext=(2,0.4),fontsize=18)
plt.annotate(r'$\tau_s = '+str(ts_s)+'s$', xy=(2,0.2),xytext=(2,0.2),fontsize=18)

ax3 = fig.add_subplot(223)
ax3.plot(tout-0.3, kkn4_d,'b-')
ax3.set_xlim([-0.3,4])
ax3.set_ylabel('Displ. or slip (m)')
ax3.set_xlabel('Time (s)')
ax3.plot(t, sd, 'k')
ax3.set_ylim([-0.05,0.6])
ax3.grid()

#Do Long Yoffe 
t,s=yoffe_reg(tr_l,ts_l)
i=where(t>=0)[0]
t=t[i]
s=s[i]
sd=cumtrapz(s,t,initial=0)
#Convolve
kkn4=convolve(s,st[0].data)*0.05
tout=arange(0,len(kkn4)*dt,dt)
kkn4_d=cumtrapz(kkn4,tout,initial=0)
kkn4=low_pass_filter(kkn4,fc,dt,2)
kkn4_d=low_pass_filter(kkn4_d,fc,dt,2)

ax2 = fig.add_subplot(222)
ax2.plot(tout-0.5, kkn4,'b-')
ax2.set_xlim([-0.3,10])
ax2.plot(t, s, 'k')
ax2.grid()
plt.annotate(r'$\tau_r = '+str(tr_l)+'s$', xy=(6,0.41),xytext=(6,0.41),fontsize=18)
plt.annotate(r'$\tau_s = '+str(ts_l)+'s$', xy=(6,0.2),xytext=(6,0.2),fontsize=18)

ax4 = fig.add_subplot(224)
ax4.plot(tout-0.5, kkn4_d,'b-')
ax4.set_xlim([-0.3,10])
ax4.set_xlabel('Time (s)')
ax4.plot(t, sd, 'k')
ax4.grid()


#ax3 = fig.add_subplot(212)
#ax3.plot(stv[0].times()-0.5, stv[0].data, 'b-')
#ax3.set_xlim([-0.5,7.5])
#ax3.set_xlabel('Time(s)')
#ax3.set_ylim([-0.3,0.4])
#ax3.set_ylabel('KKN4 vertical (m/s)')
#
#ax4 = ax3.twinx()
#ax4.plot(t_slip-0.5, slip_r, 'k')
#ax4.set_ylabel('Fault slip rate(m/s)')
#ax4.set_xlim([-0.5,7.5])
#ax4.grid()
#ax4.set_ylim([-0.9,1.2])

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(tout-0.5, kkn4,'b-')
ax1.set_xlim([-0.3,10])
ax1.set_ylabel('Velocity or slip rate (m/s)')
ax1.plot(t, s, 'k')
ax1.legend(['KKN4 vertical','Reg. Yoffe STF'])
ax1.grid()
ax1.set_xticklabels([])
plt.annotate(r'$\tau_r = '+str(tr_l)+'s$', xy=(6,0.4),xytext=(6,0.4),fontsize=18)
plt.annotate(r'$\tau_s = '+str(ts_l)+'s$', xy=(6,0.2),xytext=(6,0.2),fontsize=18)

ax3 = fig.add_subplot(212)
ax3.plot(tout-0.3, kkn4_d,'b-')
ax3.set_xlim([-0.3,10])
ax3.set_ylabel('Displacement or slip (m)')
ax3.set_xlabel('Time (s)')
ax3.plot(t, sd, 'k')
ax3.grid()

plt.subplots_adjust(left=0.2, bottom=0.15, right=0.9, top=0.98, wspace=0, hspace=0.09)

plt.show()



