from numpy.random import normal
from mtspec import mtspec


#####Run time parameters

duration=48*3600 #in seconds
dt=60.
fcorner=[1./(2*3600),1./(10*60)]   #defines the tsunami band



#####





def apply_spectrum(gauge,dt,E0=6e-7):
    '''
    Apply the modeled spectrum to the windowed time series
    '''
    
    from numpy import fft,angle,cos,sin,sqrt,mean,zeros
    from scipy.interpolate import interp1d
    
    #to frequency domain
    fourier=fft.fft(gauge)
    freq=fft.fftfreq(len(gauge),dt)
    
    #Get positive frequencies
    Nf=len(freq)
    positive_freq=freq[1:Nf/2]
    negative_freq=freq[Nf/2:]
        
    #Make POWER spectrum of white noise time series have a mean of 1
    norm_factor=mean(abs(fourier)**2)**0.5
    fourier=fourier/norm_factor
    
    #Keep phase
    phase=angle(fourier)
    
    
    #Generate long period decay
    kappa=0.4
    #P=exp(pi*kappa*f)
    
    #Generate model amplitude spectrum
    S0_positive=E0/positive_freq**2
    S0_negative=E0/negative_freq**2
    
    #Place in correct order A[0] is DC value then icnreasing positive freq then decreasing negative freq
    amplitude=zeros(len(freq))
    #DC value
    amplitude[0]=0
    #Positive freq. div by 2 to keep right power
    amplitude[1:Nf/2]=S0_positive/2
    #Negative freq
    amplitude[Nf/2:]=S0_negative/2
    
    #Apply model amplitude spectrum
    amplitude=amplitude*abs(fourier)
    
    #Obtain complex foureier series
    R=amplitude*cos(phase)
    I=amplitude*sin(phase)
    fourier=R+I*1j
    
    #ifft
    ocean_noise=fft.ifft(fourier)
    ocean_noise=ocean_noise*len(ocean_noise)
    
    return ocean_noise  




def bandpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'bandpass')
    data_filt=filtfilt(b,a,data)
    return data_filt




Nsamples=int(duration/dt)
gauge=normal(0,1,Nsamples)

#Apply background spectrum
ocean_noise=apply_spectrum(gauge,dt)

#Filter to tsunami band
ocean_noise=bandpass(ocean_noise,fcorner,1./dt,2)
psd, f, jackknife, _, _ = mtspec(data=ocean_noise, delta=dt, time_bandwidth=5,number_of_tapers=8, nfft=Nsamples, statistics=True)