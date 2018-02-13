from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime
from datetime import timedelta
from numpy import diff
from mudpy.forward import lowpass

time_epi=UTCDateTime('2015-04-25T06:11:26')


tr=10
dvert=0.5
u1=read(u'/Users/dmelgar/Slip_inv/Nepal_ttests_'+str(tr)+'/output/forward_models/'+str(tr)+'s_vr2.8.KKN4.vel.u')
u2=read(u'/Users/dmelgar/Slip_inv/Nepal_ttests_'+str(tr)+'/output/forward_models/'+str(tr)+'s_vr3.0.KKN4.vel.u')
u3=read(u'/Users/dmelgar/Slip_inv/Nepal_ttests_'+str(tr)+'/output/forward_models/'+str(tr)+'s_vr3.2.KKN4.vel.u')
u4=read(u'/Users/dmelgar/Slip_inv/Nepal_ttests_'+str(tr)+'/output/forward_models/'+str(tr)+'s_vr3.4.KKN4.vel.u')
u5=read(u'/Users/dmelgar/Slip_inv/Nepal_ttests_'+str(tr)+'/output/forward_models/'+str(tr)+'s_vr3.6.KKN4.vel.u')
u=read( u'/Users/dmelgar/Nepal2015/GPS/PPP/KKN4.LXZ.sac')

#trim
delay=15
t1=time_epi+timedelta(seconds=delay)
t2=t1+timedelta(seconds=70)
u1[0].trim(starttime=t1,endtime=t2)
u2[0].trim(starttime=t1,endtime=t2)
u3[0].trim(starttime=t1,endtime=t2)
u4[0].trim(starttime=t1,endtime=t2)
u5[0].trim(starttime=t1,endtime=t2)
u[0].trim(starttime=t1+timedelta(seconds=1),endtime=t2)
u[0].data=u[0].data-u[0].data[0]
u[0].data=lowpass(u[0].data,1.0,5.0,2)
u[0].data=diff(u[0].data/0.2)

plt.figure(figsize=(6,5))
plt.plot(u1[0].times()+delay,u1[0].data+dvert)
plt.plot(u2[0].times()+delay,u2[0].data+2*dvert)
plt.plot(u3[0].times()+delay,u3[0].data+3*dvert)
plt.plot(u4[0].times()+delay,u4[0].data+4*dvert)
plt.plot(u5[0].times()+delay,u5[0].data+5*dvert)
plt.plot(u[0].times()+delay,u[0].data,'k',lw=1.5)
plt.xlim([15,80])
plt.ylim([-0.2,3.5])
plt.grid()
plt.xlabel('Seconds after OT')
plt.ylabel('Vertical velocity (m/s)')
plt.title(str(tr)+'s rise time')

plt.annotate('Observed',xy=(65,0.1))
plt.annotate(r'$v_r=2.8km/s$',xy=(62,0.1+dvert))
plt.annotate(r'$v_r=3.0km/s$',xy=(62,0.1+2*dvert))
plt.annotate(r'$v_r=3.2km/s$',xy=(62,0.1+3*dvert))
plt.annotate(r'$v_r=3.4km/s$',xy=(62,0.1+4*dvert))
plt.annotate(r'$v_r=3.6km/s$',xy=(62,0.1+5*dvert))

plt.show()


