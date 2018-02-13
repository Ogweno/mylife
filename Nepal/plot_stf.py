from numpy import arange,where,ones,r_,diff
import matplotlib.pyplot as plt
from obspy import read

#Displacement
rise_time=1.0
rise_delta=0.5
t_slip=arange(0,10,0.01)
slip=ones(len(t_slip))
i=where(t_slip<rise_delta)[0]
slip[i]=0
i=where((t_slip>=rise_delta) & (t_slip<=(rise_delta+rise_time)))[0]
slip[i]=(1./rise_time)*t_slip[i]-rise_delta
slip_r=r_[0,diff(slip)/(t_slip[1]-t_slip[0])]

def low_pass_filter(tr,fcorner,dt,order):
    from scipy.signal import filtfilt,butter
    fnyquist=1./(2*dt)
    print fnyquist
    Fc=fcorner/fnyquist
    print Fc
    b, a = butter(order, Fc,btype='lowpass')
    y = filtfilt(b, a, tr)
    return y

st=read(u'/Users/dmelgar/Slip_inv/Nepal_Avouac_STF/GFs/dynamic/avouac.mod_13.0920.sub0011/KKN4.subfault0011.DS.disp.z')
stv=read(u'/Users/dmelgar/Slip_inv/Nepal_Avouac_STF/GFs/dynamic/avouac.mod_13.0920.sub0011/KKN4.subfault0011.DS.vel.z')
#stv[0].data=low_pass_filter(stv[0].data,5,stv[0].stats.delta,2)


fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(st[0].times()-0.5, st[0].data, 'b-')
ax1.set_xlim([-0.5,7.5])
ax1.set_ylabel('KKN4 vertical (m)')

ax2 = ax1.twinx()
ax2.plot(t_slip-0.5, slip, 'k')
ax2.set_ylabel('Fault slip (m)')
ax2.set_ylim([-0.5,2])
ax2.set_xlim([-0.5,7.5])
ax2.grid()



ax3 = fig.add_subplot(212)
ax3.plot(stv[0].times()-0.5, stv[0].data, 'b-')
ax3.set_xlim([-0.5,7.5])
ax3.set_xlabel('Time(s)')
ax3.set_ylim([-0.3,0.4])
ax3.set_ylabel('KKN4 vertical (m/s)')

ax4 = ax3.twinx()
ax4.plot(t_slip-0.5, slip_r, 'k')
ax4.set_ylabel('Fault slip rate(m/s)')
ax4.set_xlim([-0.5,7.5])
ax4.grid()
ax4.set_ylim([-0.9,1.2])
plt.show()
