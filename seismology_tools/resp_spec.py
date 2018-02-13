'''
PyResp

D. Melgar
UC Berkeley
05/2015

Compute response spectra for any kind of forcing (displacememnt, velocity or 
acceleration) and output spectral displacement, velocity, or acceleration.
This code uses second order accurate finite differences.

09/2016 Note the long period PSA has a spurios offset, I'm hunting for a fix. 
Meanwhile , I added functions from arkotke's nice pyrotd as a temporary solution
https://github.com/arkottke/pyrotd
'''

def resp_spec(forcing,frequency_vector,forcing_type='d',output_type='d',damping=5):
    '''
    forcing - ObsPy stream object with forcing time series
    frequency_vector - NumPy array with frequencies at which output is desired
    forcing_type - Type of forcing time series, either 'd','v' or 'a'
    output_type - Type of output desired, either 'd','v' or 'a'
    damping - System damping in percent.
    
    returns response,time_history     
        
    response spectra values and time histories of the oscillator at all input 
    frequencies
    '''
    
    from numpy import zeros,r_,pi,arange,ones,diff,size
    
    #Extract data
    frc=forcing[0].data
    #Sampling interval
    dt=forcing[0].stats.delta
    #How many data points?
    N=forcing[0].stats.npts
    #pad for long period spectra
    if forcing_type.lower()=='d':
        frc=r_[zeros(N),frc,frc[-1]*ones(N)]
    else:
        frc=r_[zeros(N),frc,zeros(N)]
    #Retain which indices are the non-zero-padded data
    i=arange(N,N*2)
    #How many frequencies?
    Nf=size(frequency_vector)
    #Convert damping to float value
    damping=damping*0.01
    #Set up the run
    response=zeros(Nf)
    omega=2*pi*frequency_vector #omega. Here w has Nf values
    time_history=zeros((Nf,N))
    u=zeros(N*3)
    for k in range(Nf):
        omega_current=omega[k]
        #Assign some constants
        c=(1./dt**2)+(damping*omega_current)/dt
        #c=(omega_current+((1/dt)*omega_current*damping)-(1/dt**2))
        #Solve differential equation
        if forcing_type.lower()=='a':
            for kaccel in range(2,N*3):
                u[kaccel]=(1/c)*((2*u[kaccel-1]-u[kaccel-2])/(dt**2)+((damping*omega_current*u[kaccel-2])/dt)
                    -u[kaccel-1]*omega_current**2-frc[kaccel-1])
                #u[kaccel]=(1/c)*(-frc[kaccel]+(omega_current*damping/dt)*u[kaccel-1]+
                #    (1/dt**2)*(u[kaccel-2]-2*u[kaccel-1]))
        if forcing_type.lower()=='v':
            pass
        if forcing_type.lower()=='d':
            for kdisp in range(2,N*3):
               u[kdisp]=(1/c)*(((2*u[kdisp-1]-u[kdisp-2])/(dt**2)+((damping*omega_current*u[kdisp-2])/dt)
                -u[kdisp-1]*omega_current**2+(-frc[kdisp]+2*frc[kdisp-1]-frc[kdisp-2])/(dt**2)))
               #u[kdisp]=(1/c)*((frc[kdisp-2]-2*frc[kdisp-1]+frc[kdisp])/(dt**2)+(omega_current*damping/dt)*u[kdisp-1]+
               #     (1/dt**2)*(u[kdisp-2]-2*u[kdisp-1]))
        if output_type.lower()=='v':
            u=r_[0,diff(u)]/dt
        if output_type.lower()=='a':
            u=r_[0,0,diff(u,n=2)]/(dt**2)
        #Save for output
        time_history[k,:]=u[i]
        response[k]=abs(u[i]).max()
    return response,time_history
 
    
#Stuff from pyrotd      
             
def oscillatorTimeSeries(freq, fourierAmp, oscFreq, oscDamping):
    '''Compute the time series response of an oscillator.

    Parameters
    ----------
    freq: numpy.array
        frequency of the Fourier acceleration spectrum [Hz]
    fourierAmp: numpy.array
        Fourier acceleration spectrum [g-sec]
    oscFreq: float
        frequency of the oscillator [Hz]
    oscDamping: float
        damping of the oscillator [decimal]

    Returns
    -------
    response: numpy.array
        time series response of the oscillator
    '''
    
    from numpy import power,fft
    
    # Single-degree of freedom transfer function
    h = (-power(oscFreq, 2.)
            / ((power(freq, 2.) - power(oscFreq, 2.))
                - 2.j * oscDamping * oscFreq * freq))

    # Adjust the maximum frequency considered. The maximum frequency is 5
    # times the oscillator frequency. This provides that at the oscillator
    # frequency there are at least tenth samples per wavelength.
    n = len(fourierAmp)
    m = max(n, int(2. * oscFreq / freq[1]))
    scale = float(m) / float(n)

    # Scale factor is applied to correct the amplitude of the motion for the
    # change in number of points
    return scale * fft.irfft(fourierAmp * h, 2 * (m-1))


def peakResponse(resp):
    '''Compute the maximum absolute value of a response.

    Parameters
    ----------
    resp: numpy.array
        time series of a response

    Returns
    -------
    peakResponse: float
        peak response
    '''
    from numpy import max,abs
    
    return max(abs(resp))


def responseSpectrum(timeStep, accelTs, oscFreqs, oscDamping=0.05):
    '''Compute the response spectrum for a time series.

    Parameters
    ----------
    timeStep: float
        time step of the time series [s]
    accelTs: numpy.array
        acceleration time series [g]
    oscFreqs: numpy.array
        natural frequency of the oscillators [Hz]
    oscDamping: float
        damping of the oscillator [decimal]. Default of 0.05 (i.e., 5%)

    Returns
    -------
    oscResp: numpy.array
        computed psuedo-spectral acceleartion [g]
    '''
    
    from numpy import linspace,fft,array
    
    fourierAmp = fft.rfft(accelTs)
    freq = linspace(0, 1./(2 * timeStep), num=fourierAmp.size)

    psa = [peakResponse(oscillatorTimeSeries(freq, fourierAmp, of, oscDamping))
            for of in oscFreqs]

    return array(psa)
    

