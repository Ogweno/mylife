from obspy import read
from numpy import r_,ones

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

#sta='cn23'
#delta=-10.
#e=read(path+sta+'.disp.e')
#n=read(path+sta+'.disp.n')
#u=read(path+sta+'.disp.u')
#e[0].data=delay_st(e,delta)
#n[0].data=delay_st(n,delta)
#u[0].data=delay_st(u,delta)
#e.write(outpath+sta+'.disp.e',format='SAC')
#n.write(outpath+sta+'.disp.n',format='SAC')
#u.write(outpath+sta+'.disp.u',format='SAC')
#
#
#sta='cn25'
#delta=-3.
#e=read(path+sta+'.disp.e')
#n=read(path+sta+'.disp.n')
#u=read(path+sta+'.disp.u')
#e[0].data=delay_st(e,delta)
#n[0].data=delay_st(n,delta)
#u[0].data=delay_st(u,delta)
#e.write(outpath+sta+'.disp.e',format='SAC')
#n.write(outpath+sta+'.disp.n',format='SAC')
#u.write(outpath+sta+'.disp.u',format='SAC')
#
#sta='uxal'
#delta=-8.
#e=read(path+sta+'.disp.e')
#n=read(path+sta+'.disp.n')
#u=read(path+sta+'.disp.u')
#e[0].data=delay_st(e,delta)
#n[0].data=delay_st(n,delta)
#u[0].data=delay_st(u,delta)
#e.write(outpath+sta+'.disp.e',format='SAC')
#n.write(outpath+sta+'.disp.n',format='SAC')
#u.write(outpath+sta+'.disp.u',format='SAC')
#
#sta='tnmq'
#delta=-1.
#e=read(path+sta+'.disp.e')
#n=read(path+sta+'.disp.n')
#u=read(path+sta+'.disp.u')
#e[0].data=delay_st(e,delta)
#n[0].data=delay_st(n,delta)
#u[0].data=delay_st(u,delta)
#e.write(outpath+sta+'.disp.e',format='SAC')
#n.write(outpath+sta+'.disp.n',format='SAC')
#u.write(outpath+sta+'.disp.u',format='SAC')



sta='PANG'
delta=-1.5
e=read(path+sta+'.vel.e')
n=read(path+sta+'.vel.n')
u=read(path+sta+'.vel.u')
e[0].data=delay_st(e,delta)
n[0].data=delay_st(n,delta)
u[0].data=delay_st(u,delta)
e.write(outpath+sta+'.vel.e',format='SAC')
n.write(outpath+sta+'.vel.n',format='SAC')
u.write(outpath+sta+'.vel.u',format='SAC')

sta='CMIG'
delta=3.0
e=read(path+sta+'.vel.e')
n=read(path+sta+'.vel.n')
u=read(path+sta+'.vel.u')
e[0].data=delay_st(e,delta)
n[0].data=delay_st(n,delta)
u[0].data=delay_st(u,delta)
e.write(outpath+sta+'.vel.e',format='SAC')
n.write(outpath+sta+'.vel.n',format='SAC')
u.write(outpath+sta+'.vel.u',format='SAC')

sta='TGIG'
delta=1.5
e=read(path+sta+'.vel.e')
n=read(path+sta+'.vel.n')
u=read(path+sta+'.vel.u')
e[0].data=delay_st(e,delta)
n[0].data=delay_st(n,delta)
u[0].data=delay_st(u,delta)
e.write(outpath+sta+'.vel.e',format='SAC')
n.write(outpath+sta+'.vel.n',format='SAC')
u.write(outpath+sta+'.vel.u',format='SAC')


sta='TNMQ'
delta=-3.5
e=read(path+sta+'.disp.e')
n=read(path+sta+'.disp.n')
u=read(path+sta+'.disp.u')
e[0].data=delay_st(e,delta)
n[0].data=delay_st(n,delta)
u[0].data=delay_st(u,delta)
e.write(outpath+sta+'.disp.e',format='SAC')
n.write(outpath+sta+'.disp.n',format='SAC')
u.write(outpath+sta+'.disp.u',format='SAC')


sta='OXUM'
delta=-2.0
e=read(path+sta+'.disp.e')
n=read(path+sta+'.disp.n')
u=read(path+sta+'.disp.u')
e[0].data=delay_st(e,delta)
n[0].data=delay_st(n,delta)
u[0].data=delay_st(u,delta)
e.write(outpath+sta+'.disp.e',format='SAC')
n.write(outpath+sta+'.disp.n',format='SAC')
u.write(outpath+sta+'.disp.u',format='SAC')


sta='TGIG'
delta=-2.0
e=read(path+sta+'.disp.e')
n=read(path+sta+'.disp.n')
u=read(path+sta+'.disp.u')
e[0].data=delay_st(e,delta)
n[0].data=delay_st(n,delta)
u[0].data=delay_st(u,delta)
e.write(outpath+sta+'.disp.e',format='SAC')
n.write(outpath+sta+'.disp.n',format='SAC')
u.write(outpath+sta+'.disp.u',format='SAC')