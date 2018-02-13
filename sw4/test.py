from obspy import read
from scipy.integrate import cumtrapz
from matplotlib import pyplot as plt
from numpy import r_,diff


def lowpass(data, fcorner, fsample, order, zerophase=True):
    '''
    Make a highpass zero phase filter
    '''
    from scipy.signal import butter, filtfilt, lfilter
    from numpy import size, array

    fnyquist = fsample / 2
    b, a = butter(order, array(fcorner) / (fnyquist), 'lowpass')
    if zerophase == True:
        data_filt = filtfilt(b, a, data)
    else:
        data_filt = lfilter(b, a, data)

    return data_filt




d=read(u'/Users/dmelgar/code/sw4/examples/rupture/napa_ramp/S3.y')
v=read(u'/Users/dmelgar/code/sw4/examples/rupture/napa_triangle/S3.y')

ddiff=r_[0,diff(d[0].data)/d[0].stats.delta]
vfilt=lowpass(v[0].data, 0.5, 1./v[0].stats.delta, 2)
vint=cumtrapz(v[0].data,v[0].times(),initial=0)


plt.figure()
plt.plot(d[0].times(),d[0].data)
plt.plot(v[0].times(),v[0].data)
plt.legend(['ramp','triangle'])


plt.figure()
plt.plot(d[0].times(),d[0].data)
plt.plot(v[0].times(),vint)
plt.legend(['ramp','integrated triangle'])


plt.figure()
plt.plot(v[0].times(),v[0].data)
plt.plot(d[0].times(),ddiff)
plt.legend(['triangle','differentiated ramp'])

plt.figure()
plt.plot(v[0].times(),v[0].data)
plt.plot(v[0].times(),vfilt)
plt.legend(['velocity','filtered velocity'])

plt.show()