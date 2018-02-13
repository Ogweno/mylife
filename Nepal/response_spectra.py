from obspy import read
from resp_spec import resp_spec
from numpy import interp,arange,logspace,log10,c_,savetxt

dt=0.01
frequency_vector=logspace(log10(1./50),log10(20),100)

sta='NAST'
e=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXE.sac')
n=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXN.sac')
u=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXZ.sac')
#Resample to 100Hz
t=arange(0,e[0].times().max()+dt,dt)
e[0].data=interp(t,e[0].times(),e[0].data)
e[0].stats.delta=dt
n[0].data=interp(t,n[0].times(),n[0].data)
n[0].stats.delta=dt
u[0].data=interp(t,u[0].times(),u[0].data)
u[0].stats.delta=dt
#Response
Re,time_history=resp_spec(e,frequency_vector,forcing_type='d',output_type='a',damping=5.)
Rn,time_history=resp_spec(n,frequency_vector,forcing_type='d',output_type='a',damping=5.)
Ru,time_history=resp_spec(e,frequency_vector,forcing_type='d',output_type='a',damping=5.)
savetxt(u'/Users/dmelgar/Nepal2015/Response_Spectra/'+sta+'.acc.resp',c_[frequency_vector,Rn,Re,Ru],fmt='%.6f\t%.6f\t%.6f\t%.6f')

sta='KKN4'
e=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXE.sac')
n=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXN.sac')
u=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXZ.sac')
#Resample to 100Hz
t=arange(0,e[0].times().max()+dt,dt)
e[0].data=interp(t,e[0].times(),e[0].data)
e[0].stats.delta=dt
n[0].data=interp(t,n[0].times(),n[0].data)
n[0].stats.delta=dt
u[0].data=interp(t,u[0].times(),u[0].data)
u[0].stats.delta=dt
#Response
Re,time_history=resp_spec(e,frequency_vector,forcing_type='d',output_type='a',damping=5.)
Rn,time_history=resp_spec(n,frequency_vector,forcing_type='d',output_type='a',damping=5.)
Ru,time_history=resp_spec(e,frequency_vector,forcing_type='d',output_type='a',damping=5.)
savetxt(u'/Users/dmelgar/Nepal2015/Response_Spectra/'+sta+'.acc.resp',c_[frequency_vector,Rn,Re,Ru],fmt='%.6f\t%.6f\t%.6f\t%.6f')

sta='KATNP'
e=read(u'/Users/dmelgar/Nepal2015/strong_motion/'+sta+'.acc.e')
n=read(u'/Users/dmelgar/Nepal2015/strong_motion/'+sta+'.acc.n')
u=read(u'/Users/dmelgar/Nepal2015/strong_motion/'+sta+'.acc.u')
#Response
Re,time_history=resp_spec(e,frequency_vector,forcing_type='a',output_type='a',damping=5.)
Rn,time_history=resp_spec(n,frequency_vector,forcing_type='a',output_type='a',damping=5.)
Ru,time_history=resp_spec(e,frequency_vector,forcing_type='a',output_type='a',damping=5.)
savetxt(u'/Users/dmelgar/Nepal2015/Response_Spectra/'+sta+'.acc.resp',c_[frequency_vector,Rn,Re,Ru],fmt='%.6f\t%.6f\t%.6f\t%.6f')
