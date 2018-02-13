import resp_spec
from obspy import read,Stream,Trace
from numpy import logspace
from scipy.integrate import cumtrapz
from matplotlib import pyplot as plt

#Read  Data
a=read(u'/Users/dmelgar/Amatrice2016/strong_motion/sac/NRC.HNN.sac')

#Trim cause I'm lazy
a[0].trim(endtime=a[0].stats.starttime+40)

#Freq for Sa
f=logspace(-2,1,100)

#Filter corner for itnegrating strong motion
fcorner=1./0.7

#Define a high pass filter
def highpass_filter(data,fcorner,fsample,order):
    '''
    Make a highpass zero phase filter
    '''
    from numpy import size,array
    from scipy.signal import butter,filtfilt
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt

#De-mean and high pass filter accel
a[0].data=a[0].data-a[0].data.mean()
a[0].data=highpass_filter(a[0].data,fcorner,1./a[0].stats.delta,2)

#Integrate to displacememnt
d=Stream(Trace())
v=cumtrapz(a[0].data,a[0].times(),initial=0)
d[0].data=cumtrapz(v,a[0].times(),initial=0)
d[0].stats.delta=a[0].stats.delta

#Plot to make sure no funny bussiness
plt.figure()
plt.subplot(211)
plt.plot(a[0].times(),a[0].data,'k')
plt.legend(['Acc (m/s/s)'])

plt.subplot(212)
plt.plot(d[0].times(),d[0].data,'r')
plt.xlabel('Seconds')
plt.legend(['Disp (m)'])


#Run resp_spec on acceleration and displ and collect Sa output
Sa_from_accel,thist=resp_spec.resp_spec(a,f,forcing_type='a',output_type='a',damping=5.)
Sa_from_disp,thist=resp_spec.resp_spec(d,f,forcing_type='d',output_type='a',damping=5.)

plt.figure()
plt.plot(1./f,Sa_from_accel,'k')
plt.plot(1./f,Sa_from_disp,'r')
plt.xlim([0,10])
plt.legend(['From accel','From disp'])
plt.xlabel('Period (s)')
plt.ylabel('Sa (m/s/s)')

#Sd_from_accel,thist=resp_spec.resp_spec(a,f,forcing_type='a',output_type='d',damping=5.)
#Sd_from_disp,thist=resp_spec.resp_spec(d,f,forcing_type='d',output_type='d',damping=5.)
#
#plt.figure()
#plt.plot(1./f,Sd_from_accel,'k')
#plt.plot(1./f,Sd_from_disp,'r')
#plt.xlim([0,10])
#plt.legend(['From accel','From disp'])
#plt.xlabel('Period (s)')
#plt.ylabel('Sd (m)')



#test pyrotd stuff

Sa_pyrotd=resp_spec.responseSpectrum(a[0].stats.delta, a[0].data, f)

plt.figure()
plt.plot(1./f,Sa_from_accel,'k')
plt.plot(1./f,Sa_pyrotd,'r')
plt.xlim([0,10])
plt.legend(['Diego''s crap code','pyrotd'])
plt.xlabel('Period (s)')
plt.ylabel('Sa (m/s/s)')

plt.show()