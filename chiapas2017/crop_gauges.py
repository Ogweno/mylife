from obspy import read
from numpy import array
from datetime import timedelta

sta=['43413','sali','huat','pchi','ptan']
tcrop=array([10000,8944,5495,9400,5879])

for k in range(len(sta)):
    st=read(u'/Users/dmelgar/Chiapas2017/tsunami/4inversion/'+sta[k]+'.sac')
    st[0].trim(endtime=st[0].stats.starttime+timedelta(seconds=tcrop[k]))
    st.write(u'/Users/dmelgar/Chiapas2017/tsunami/4inversion/trim/'+sta[k]+'.sac',format='SAC')