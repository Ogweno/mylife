from obspy import read
from numpy import r_,ones,zeros

path=u'/Users/dmelgar/Slip_inv/Chiapas_hernandez_new/data/waveforms/before_delta_t/'
outpath='/Users/dmelgar/Slip_inv/Chiapas_hernandez_new/data/waveforms/'


def delay_st(st,delta):
    d=st[0].data
    npts=int(abs(delta)/st[0].stats.delta)
    if delta<0:
        d=r_[d[npts:-1],d[-1]*ones(npts+1)]
    else:
        d=r_[ones(npts)*d[0],d[0:-npts]]
    return d

sta='43413'
delta=-3*60.
e=read(path+sta+'.sac')
e[0].data=delay_st(e,delta)
e.write(outpath+sta+'.sac',format='SAC')

sta='huat'
delta=-2*60.
e=read(path+sta+'.sac')
e[0].data=delay_st(e,delta)
e.write(outpath+sta+'.sac',format='SAC')


sta='ptan'
delta=-2*60.
e=read(path+sta+'.sac')
e[0].data=delay_st(e,delta)
e.write(outpath+sta+'.sac',format='SAC')


