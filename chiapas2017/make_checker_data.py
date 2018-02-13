from numpy import genfromtxt,load,zeros,arange,save,r_,array,squeeze
from mudpy.inverse import ds2rot

f=genfromtxt('/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/forward_models/checker_new.rupt')

#fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/static_1win_data.npy"
#G=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/static_1win.npy")

#fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/hrgps_1win_data.npy"
#G=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/hrgps_1win.npy")

#fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/sm_1win_data.npy"
#G=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/sm_1win.npy")

fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/tsunami_1win_data.npy"
G=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/tsunami_1win.npy")

#fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/static_hrgps_sm_tsunami_1win_data.npy"
#foutG="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/static_hrgps_sm_tsunami_1win.npy"
#G1=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/static_1win.npy")
#G2=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/hrgps_1win.npy")
#G3=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/sm_1win.npy")
#G4=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/tsunami_1win.npy")
#G=r_[G1,G2,G3,G4]
#save(foutG,G)


#fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/gps_onewindow_data.npy"
#G=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/gps_1win_noshallow_vr3.6.npy")

##
#fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/tsunami_onewindow_data.npy"
#G=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/tsunami_1win_noshallow_vr3.6.npy")

#fout="/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/static_gps_tsunami_onewindow_data.npy"
#G1=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/static_1win_noshallow_vr3.6.npy")
#G2=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/gps_1win_noshallow_vr3.6.npy")
#G3=load("/Users/dmelgar/Slip_inv/Chiapas_hernandez_checker/GFs/matrices/tsunami_1win_noshallow_vr3.6.npy")
#G=r_[G1,G2,G3]



def lowpass(data,fcorner,fsample,order,zerophase=True):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt,lfilter
    from numpy import size,array
    
    if size(fcorner)==2:
        ftype='bandpass'
    else:
        ftype='lowpass'
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),ftype)
    if zerophase==True:
        data_filt=filtfilt(b,a,data)
    else:
        data_filt=lfilter(b,a,data)
        
    return data_filt






m=zeros((len(f)*2,1))
i=arange(1,len(f)*2,2)
m[i,0]=f[:,9]
mrot=ds2rot(m,225)

d=G.dot(mrot)


#sm
#d=lowpass(squeeze(d),array([0.02,0.1]),4,2)

#tsunami
#d=lowpass(squeeze(d),array([1./7200,1./600]),1./60,2)

save(fout,d)
