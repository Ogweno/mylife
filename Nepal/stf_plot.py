from numpy import arange,pi,where,convolve,r_,diff,zeros
from matplotlib import pyplot as plt
from obspy import read
import matplotlib

tr=3.3 #How long it lasts
ts=1.7 #How fast it rises
kkn4_offset=21.5
katnp_offset=21.5
nast_offset=23.3

katnp_dy=0.2
nast_dy=0.12

dt=0.01
t=arange(dt,100,dt)
tstart=4
max_sr=0.58
matplotlib.rcParams.update({'font.size': 16})

#Load up some data
katnp=read(u'/Users/dmelgar/Slip_inv/Nepal_Avouac/data/waveforms/KATNP.vel.u')
nast=read(u'/Users/dmelgar/Slip_inv/Nepal_Avouac/data/waveforms/NAST.disp.u')
kkn4=read(u'/Users/dmelgar/Slip_inv/Nepal_Avouac/data/waveforms/KKN4.disp.u')
nast[0].data=r_[0,diff(nast[0].data)/nast[0].stats.delta]
kkn4[0].data=r_[0,diff(kkn4[0].data)/kkn4[0].stats.delta]

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
s=(s[0:len(t)]/s.max())*max_sr
s=r_[zeros(60/dt),s]
t=arange(-60,len(s)*dt-60,dt)
#plt.figure()
#plt.plot(s)
#plt.show()

plt.figure()
plt.plot(kkn4[0].times(),kkn4[0].data+1.2,'#606060')
plt.plot(katnp[0].times(),katnp[0].data+0.6,'#DAA520')
plt.plot(nast[0].times(),nast[0].data,'#DC143C')
plt.plot(t+kkn4_offset,s+1.2,'#1E90FF')aplt.plot(t+nast_offset,s-nast_dy,'#1E90FF')
plt.legend(['KKN4','KATNP','NAST','Tinti STF'],bbox_to_anchor=(0.9,1.1),ncol=4,frameon=False,fontsize=15)
plt.ylabel('Up(m/s)')
plt.xlabel('Seconds after OT')
plt.xlim([10,60])
plt.grid()
plt.annotate(r'$\tau_r = 3.3s$', xy=(45,1.7),xytext=(45,1.7),fontsize=18)
plt.annotate(r'$\tau_s = 1.7s$', xy=(45,1.55),xytext=(45,1.55),fontsize=18)
plt.show()
