'''
D.Melgar
Decimate Nepal waveforms
'''

from glob import glob
from obspy import read
from scipy.signal import butter,filtfilt

def stdecimate(st,factor,order):
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout

#list files to decimate
files=glob('/Users/dmelgar/Slip_inv/nepal_example_ch/data/waveforms/*')
for k in range(len(files)):
    print files[k]
    st=read(files[k])
    #decimate it
    stdec=stdecimate(st,5,4)
    #save
    stdec.write(files[k],format='SAC')
    