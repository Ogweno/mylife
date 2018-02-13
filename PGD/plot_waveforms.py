from matplotlib import pyplot as plt
from obspy import read
from numpy import mean
from obspy.core import UTCDateTime
from datetime import timedelta

def highpass(data,fcorner,fsample,order=5):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt

f1=u'/Users/dmelgar/Maule2010/GPS/proc/cons.LXE.sac'
f2=u'/Users/dmelgar/Iquique2014/GPS/proc/psga.LXE.sac'
f3=u'/Users/dmelgar/Nicoya2012/GPS/proc/caba.LXZ.sac'
f4=u'/Users/dmelgar/Napa2014/GPS/sac/p267.LXN.sac'
t1=UTCDateTime('2010-02-27T06:34:14')
t2=UTCDateTime('2014-04-01T23:46:47')
t3=UTCDateTime('2012-09-05T14:42:08')
t4=UTCDateTime('2014-08-24T10:20:44')

w1=read(f1)
w1[0].data=w1[0].data-mean(w1[0].data[23000:23500])
w11=w1.copy()
w12=w1.copy()
w13=w1.copy()
w14=w1.copy()

w2=read(f2)
w2[0].data=w2[0].data-mean(w2[0].data[85000:85500])
w21=w2.copy()
w22=w2.copy()
w23=w2.copy()
w24=w2.copy()

w3=read(f3)
w3[0].data=w3[0].data-mean(w3[0].data[241945:263485])
w31=w3.copy()
w32=w3.copy()
w33=w3.copy()
w34=w3.copy()

w4=read(f4)
w4[0].data=w4[0].data-mean(w4[0].data[500:600])
w41=w4.copy()
w42=w4.copy()
w43=w4.copy()
w44=w4.copy()

#Filter
w11[0].data=highpass(w1[0].data,0.02,1./w1[0].stats.delta)
w12[0].data=highpass(w1[0].data,0.07,1./w1[0].stats.delta)
w13[0].data=highpass(w1[0].data,0.1,1./w1[0].stats.delta)
w14[0].data=highpass(w1[0].data,0.2,1./w1[0].stats.delta)

w21[0].data=highpass(w2[0].data,0.02,1./w2[0].stats.delta)
w22[0].data=highpass(w2[0].data,0.07,1./w2[0].stats.delta)
w23[0].data=highpass(w2[0].data,0.1,1./w2[0].stats.delta)
w24[0].data=highpass(w2[0].data,0.2,1./w2[0].stats.delta)

w31[0].data=highpass(w3[0].data,0.02,1./w3[0].stats.delta)
w32[0].data=highpass(w3[0].data,0.07,1./w3[0].stats.delta)
w33[0].data=highpass(w3[0].data,0.1,1./w3[0].stats.delta)
w34[0].data=highpass(w3[0].data,0.2,1./w3[0].stats.delta)

w41[0].data=highpass(w4[0].data,0.02,1./w4[0].stats.delta)
w42[0].data=highpass(w4[0].data,0.07,1./w4[0].stats.delta)
w43[0].data=highpass(w4[0].data,0.1,1./w4[0].stats.delta)
w44[0].data=highpass(w4[0].data,0.2,1./w4[0].stats.delta)

#Trim
DT=timedelta(seconds=180)
w1.trim(starttime=t1,endtime=t1+DT)
w11.trim(starttime=t1,endtime=t1+DT)
w12.trim(starttime=t1,endtime=t1+DT)
w13.trim(starttime=t1,endtime=t1+DT)
w14.trim(starttime=t1,endtime=t1+DT)

DT=timedelta(seconds=180)
w2.trim(starttime=t2,endtime=t2+DT)
w21.trim(starttime=t2,endtime=t2+DT)
w22.trim(starttime=t2,endtime=t2+DT)
w23.trim(starttime=t2,endtime=t2+DT)
w24.trim(starttime=t2,endtime=t2+DT)

DT=timedelta(seconds=180)
w3.trim(starttime=t3,endtime=t3+DT)
w31.trim(starttime=t3,endtime=t3+DT)
w32.trim(starttime=t3,endtime=t3+DT)
w33.trim(starttime=t3,endtime=t3+DT)
w34.trim(starttime=t3,endtime=t3+DT)

DT=timedelta(seconds=180)
w4.trim(starttime=t4,endtime=t4+DT)
w41.trim(starttime=t4,endtime=t4+DT)
w42.trim(starttime=t4,endtime=t4+DT)
w43.trim(starttime=t4,endtime=t4+DT)
w44.trim(starttime=t4,endtime=t4+DT)

#Plot
plt.figure()
plt.subplot(411)
plt.plot(w1[0].times(),w1[0].data,'k',lw=2)
plt.plot(w11[0].times(),w11[0].data,'#FF6347',lw=1)
plt.plot(w12[0].times(),w12[0].data,'#1E90FF',lw=1)
plt.plot(w13[0].times(),w13[0].data,'#32CD32',lw=1)
plt.xlim([0,100])
plt.yticks([-4,-2,0])
plt.grid()
plt.ylabel('Disp. (m)')
plt.text(3,-5,'Maule, station CONS, east')
plt.text(142,-0.5,'(a)')
plt.legend(['Unfiltered',r'$f_c=0.02Hz$',r'$f_c=0.07Hz$',r'$f_c=0.1Hz$',r'$f_c=0.2Hz$'],bbox_to_anchor=(0., 1.10, 1., .102), loc=3,
                ncol=4, mode="expand", borderaxespad=0.,fontsize=12)

plt.subplot(412)
plt.plot(w2[0].times(),w2[0].data,'k',lw=2)
plt.plot(w21[0].times(),w21[0].data,'#FF6347',lw=1)
plt.plot(w22[0].times(),w22[0].data,'#1E90FF',lw=1)
plt.plot(w23[0].times(),w23[0].data,'#32CD32',lw=1)
plt.grid()
plt.ylabel('Disp. (m)')
plt.yticks([-0.8,-0.4,0])
plt.xlim([0,100])
plt.text(2,-0.8,'Iquique, station PSGA, east')
plt.text(114,-0.1,'(b)')

plt.subplot(413)
plt.plot(w3[0].times(),w3[0].data,'k',lw=2)
plt.plot(w31[0].times(),w31[0].data,'#FF6347',lw=1)
plt.plot(w32[0].times(),w32[0].data,'#1E90FF',lw=1)
plt.plot(w33[0].times(),w33[0].data,'#32CD32',lw=1)
plt.grid()
plt.ylabel('Disp. (m)')
plt.yticks([-0.2,0])
plt.text(1,0.1,'Nicoya, station CABA, vert.')
plt.text(85,0.1,'(c)')
plt.xlim([0,100])

plt.subplot(414)
plt.plot(w4[0].times(),w4[0].data,'k',lw=2)
plt.plot(w41[0].times(),w41[0].data,'#FF6347',lw=1)
plt.plot(w42[0].times(),w42[0].data,'#1E90FF',lw=1)
plt.plot(w43[0].times(),w43[0].data,'#32CD32',lw=1)
plt.xlim([0,100])
plt.yticks([-0.05,0,0.05])
plt.grid()
plt.text(1,-0.063,'Napa, station P267, north')
plt.text(85,0.05,'(d)')
plt.ylabel('Disp. (m)')
plt.xlabel('Seconds after origin time')

plt.show()