from numpy import fft,angle,cos,sin,sqrt,mean,arange,exp,pi,real

dt=0.1
t=arange(-2,2,dt)
y=exp(-pi*t**2)


#ifft test
fourier=fft.fft(y)
#fourier=abs(fourier)
f=fft.fftfreq(len(y),dt)
yout=real(fft.ifft(fourier))
yout=yout

fourier2=exp(-pi*f**2)
yout2=real(fft.ifft(fourier2))
#Ok, define 


##Make POWER spectrum of windowed time series have a mean of 1
##norm_factor=hf_dt*mean(abs(fourier)**2)**0.5
##norm_factor=mean(abs(fourier))
#norm_factor=mean(abs(fourier)**2)**0.5
#fourier=fourier/norm_factor
#
##Keep phase
#phase=angle(fourier)
#
##resample model amplitude spectr to frequencies
#interp=interp1d(f,A,bounds_error=False)
#amplitude=interp(freq)
#amplitude[0]=amplitude[1]
#
##Apply model amplitude spectrum
#amplitude=amplitude*abs(fourier)
#
##Obtain complex foureier series
#R=amplitude*cos(phase)
#I=amplitude*sin(phase)
#fourier=R+I*1j
#
##ifft
#seis=fft.irfft(fourier)
#seis=seis*len(seis)