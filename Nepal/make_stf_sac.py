from numpy import arange,zeros,pi,where,r_
from scipy.signal import convolve
from scipy.interpolate import interp1d
from obspy import Stream,Trace
from scipy.integrate import trapz
from matplotlib import pyplot as plt

fout='/Users/dmelgar/Slip_inv/Nepal_ttests_7/GFs/reg_yoffe_7.sac'
tr=4.55
ts=1.4
teff=tr+2*ts

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

def yoffe_reg(tr,ts,dt):
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

t,s=yoffe_reg(tr,ts,0.01)

#Integrate, then normalize
I=trapz(s,t)
print I
s=s/I
s=s*0.2 #Normalize to dt interval
#Cut
i=where((t>=0) & (t<(teff+1)))[0]
t=t[i]
s=s[i]
#decimate
ti=arange(0,t[-1],0.2)
f=interp1d(t,s,bounds_error=False,fill_value=0)
si=f(ti)
plt.figure()
plt.plot(t,s,ti,si)
plt.show()

#Put ins ac file
st=Stream(Trace())
st[0].data=si
st[0].stats.delta=0.2
st.write(fout,format='SAC')