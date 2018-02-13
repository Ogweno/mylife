from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from matpltolib import pyplot as plt
from obspy.geodetics import gps2dist_azimuth
'''
Afghanistan is GMT+4:30

Source is @Shigal, Afghanistan

71.229E, 35.110N

Kabul is at

69.167E 34.533N
'''

gps2dist_azimuth(35.110,71.229,34.553,69.167)
d,az,baz=gps2dist_azimuth(35.110,71.229,34.553,69.167)
p_vel=6.0 #km/s
tp=(d/1000)/p_vel
tsource=19+32/60.+tp/3600.

client = Client("IRIS")


starttime = UTCDateTime('2017-04-13T00:00:00')
endtime = UTCDateTime('2017-04-13T23:59:59.99')

st = client.get_waveforms("IU", "KBL", "", "BH?", starttime, endtime)
#hours=(st[0].times()/3600)+4.5
hours=st[0].times()

plt.figure()
plt.subplot(311)
plt.plot(hours,st[0].data)
#plt.scatter(tsource,0,marker='|',s=30,color='r')

plt.subplot(312)
plt.plot(hours,st[1].data)

plt.subplot(313)
plt.plot(hours,st[2].data)

plt.show()