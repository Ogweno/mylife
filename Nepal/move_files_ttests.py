from glob import glob
from obspy import read
files=glob(u'/Users/dmelgar/Slip_inv/Nepal_ttests_6/output/forward_models/*vr3.2*.vel*')
dirchop=u'/Users/dmelgar/Slip_inv/Nepal_Avouac_1s/data/waveforms/'
dirout='/Users/dmelgar/Slip_inv/Nepal_ttests_inv/data/waveforms/'
for k in range(len(files)):
    file_out=dirout+files[k].split('/')[-1].split('.')[2]+'.'+files[k].split('/')[-1].split('.')[3]+'.'+files[k].split('/')[-1].split('.')[4]
    file_chop=dirchop+files[k].split('/')[-1].split('.')[2]+'.'+files[k].split('/')[-1].split('.')[3]+'.'+files[k].split('/')[-1].split('.')[4]
    st=read(files[k])
    stchop=read(file_chop)
    st[0].trim(starttime=stchop[0].stats.starttime,endtime=stchop[0].stats.endtime)
    st.write(file_out,format='SAC')