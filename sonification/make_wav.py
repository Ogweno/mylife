import numpy as np
import wave
from obspy import read
from scipy.interpolate import interp1d
from scipy.io.wavfile import write


fout='GO04_mainshock_and_aftershocks_normalspeed.wav'
speedup=1

#read in raw data
#st=read(u'/Users/dmelgar/Coquimbo2015/strong_motion/proc/GO04.HNN.sac')
st=read(u'/Users/dmelgar/Coquimbo2015/strong_motion/trim/GO04.HNN.sac')

samplerate=44100

#remove baseline
st[0].data=st[0].data-np.mean(st[0].data[0:100])

#Speedup
st[0].stats.delta=st[0].stats.delta/speedup

#Resample
f=interp1d(st[0].times(),st[0].data,kind='linear')
t=np.arange(0,np.max(st[0].times()),1./samplerate)
y=f(t)

#Write
data=y
scaled = np.int16(y/np.max(np.abs(y)) * 32767)
write(fout, 44100, scaled)