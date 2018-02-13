from numpy import genfromtxt,r_,where,mean
from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime
import nitime.algorithms as tsa
from matplotlib.ticker import MultipleLocator
from mtspec import mtspec

st=read('/Users/dmelgar/tidegauge_noise/crescent_city/cres.sac')
st[0].data=st[0].data-mean(st[0].data)
Tw=6
Ntapers=20
psd, f, jackknife, _, _ = mtspec(data=st[0].data, delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)


#filter
fcorner=1./(2*3600)
def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt
st_filt=st.copy()
st_filt[0].data=highpass(st[0].data,fcorner,1./st[0].stats.delta,2)




plt.figure(figsize=(12,6))

ax=plt.axes([0.1,0.15,0.5,0.33])
ax.xaxis.set_ticks_position('both')
ax.plot(st_filt[0].times()/86400,st_filt[0].data,c='#303030',lw=0.8)
ax.set_xlim([0,30])
ax.set_ylim([-0.2,0.2])
#ax.yaxis.set_ticklabels(['','-0.1','','0','','0.1',''])
ax.set_ylabel('Sea level (m)')
ax.set_xlabel('Days')
xmajorLocator = MultipleLocator(5)
xminorLocator = MultipleLocator(1)
ymajorLocator = MultipleLocator(0.1)
yminorLocator = MultipleLocator(0.05)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.legend(['Filtered'],bbox_to_anchor=(0.7, 1), loc=2)

ax=plt.axes([0.1,0.52,0.5,0.33])
ax.xaxis.tick_top()
ax.xaxis.set_ticks_position('both')
ax.plot(st[0].times()/86400,st[0].data,c='#303030',lw=0.8)
ax.set_xlim([0,30])
ax.set_ylim([-2,2])
#ax.yaxis.set_ticklabels(['','-0.1','','0','','0.1',''])
ax.set_ylabel('Sea level (m)')
ax.set_xlabel('Days')
ax.xaxis.set_label_position("top")
xmajorLocator = MultipleLocator(5)
xminorLocator = MultipleLocator(1)
ymajorLocator = MultipleLocator(1)
yminorLocator = MultipleLocator(0.2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.legend(['Raw'],loc=0)



ax=plt.axes([0.62,0.15,0.3,0.7])
ax.yaxis.tick_right()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.plot((1./f)/60,psd,lw=0.5,c='#303030')
ax.set_xlim([12,120])
ax.set_ylim([0,12])
ax.set_ylabel(r'PSD ($m^2\times\; min$)')
ax.yaxis.set_label_position("right")
ax.set_xlabel('Period (min)')
xmajorLocator = MultipleLocator(20)
xminorLocator = MultipleLocator(5)
ymajorLocator = MultipleLocator(2)
yminorLocator = MultipleLocator(0.2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)





plt.show()