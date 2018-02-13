def bandpass(tr,fcorner,dt,order):
    from scipy.signal import filtfilt,butter
    from numpy import array
    fnyquist=1./(2*dt)
    print fnyquist
    Fc=array(fcorner)/fnyquist
    print Fc
    b, a = butter(order, Fc,btype='bandpass')
    y = filtfilt(b, a, tr)
    return y
    
def lowpass(tr,fcorner,dt,order):
    from scipy.signal import filtfilt,butter
    fnyquist=1./(2*dt)
    print fnyquist
    Fc=fcorner/fnyquist
    print Fc
    b, a = butter(order, Fc,btype='lowpass')
    y = filtfilt(b, a, tr)
    return y
    
def highpass(tr,fcorner,dt,order):
    from scipy.signal import filtfilt,butter
    fnyquist=1./(2*dt)
    print fnyquist
    Fc=fcorner/fnyquist
    print Fc
    b, a = butter(order, Fc,btype='highpass')
    y = filtfilt(b, a, tr)
    return y