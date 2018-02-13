from obspy import read
from scipy.signal import butter,filtfilt
from scipy.integrate import cumtrapz
from numpy import array
from matplotlib import pyplot as plt

n=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/PPIG.HLE.sac')
e=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/PPIG.HLN.sac')
z=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/PPIG.HLZ.sac')

n2=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/PLIG.HLE.sac')
e2=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/PLIG.HLN.sac')
z2=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/PLIG.HLZ.sac')

fcorner=0.05
fcorner2=1.0

dt1=107
xlim1=[0,90]

def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
        
    return data_filt
    
def lowpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'lowpass')
    data_filt=filtfilt(b,a,data)
        
    return data_filt
    
    
dn=n.copy()
de=e.copy()
dz=z.copy()

dn[0].data=cumtrapz(cumtrapz(n[0].data,n[0].times(),initial=0),n[0].times(),0)
de[0].data=cumtrapz(cumtrapz(e[0].data,e[0].times(),initial=0),e[0].times(),0)
dz[0].data=cumtrapz(cumtrapz(z[0].data,z[0].times(),initial=0),z[0].times(),0)

dn[0].data=highpass(dn[0].data,fcorner,1./dn[0].stats.delta,2)
de[0].data=highpass(de[0].data,fcorner,1./dn[0].stats.delta,2)
dz[0].data=highpass(dz[0].data,fcorner,1./dn[0].stats.delta,2)


lp_dn=dn.copy()
lp_de=de.copy()
lp_dz=dz.copy()

lp_dn[0].data=lowpass(dn[0].data,fcorner2,1./dn[0].stats.delta,2)
lp_de[0].data=lowpass(de[0].data,fcorner2,1./dn[0].stats.delta,2)
lp_dz[0].data=lowpass(dz[0].data,fcorner2,1./dn[0].stats.delta,2)






dn2=n2.copy()
de2=e2.copy()
dz2=z2.copy()

dn2[0].data=cumtrapz(cumtrapz(n2[0].data,n2[0].times(),initial=0),n2[0].times(),0)
de2[0].data=cumtrapz(cumtrapz(e2[0].data,e2[0].times(),initial=0),e2[0].times(),0)
dz2[0].data=cumtrapz(cumtrapz(z2[0].data,z2[0].times(),initial=0),z2[0].times(),0)

dn2[0].data=highpass(dn2[0].data,fcorner,1./dn2[0].stats.delta,2)
de2[0].data=highpass(de2[0].data,fcorner,1./dn2[0].stats.delta,2)
dz2[0].data=highpass(dz2[0].data,fcorner,1./dn2[0].stats.delta,2)


lp_dn2=dn2.copy()
lp_de2=de2.copy()
lp_dz2=dz2.copy()

lp_dn2[0].data=lowpass(dn2[0].data,fcorner2,1./dn[0].stats.delta,2)
lp_de2[0].data=lowpass(de2[0].data,fcorner2,1./dn[0].stats.delta,2)
lp_dz2[0].data=lowpass(dz2[0].data,fcorner2,1./dn[0].stats.delta,2)








plt.figure(figsize=(12,8))

plt.subplot(431)
plt.plot(n[0].times()-dt1,n[0].data*100)
plt.title('North',fontsize=14)
plt.ylabel('Accel (cm/s/s)',fontsize=14)
plt.legend(['PPIG'],fontsize=14)
plt.xlim(xlim1)
plt.subplot(432)
plt.plot(e[0].times()-dt1,e[0].data*100)
plt.title('East',fontsize=14)
plt.xlim(xlim1)
plt.subplot(433)
plt.plot(z[0].times()-dt1,z[0].data*100)
plt.title('Up',fontsize=14)
plt.xlim(xlim1)

plt.subplot(434)
plt.plot(dn[0].times()-dt1,dn[0].data*100)
plt.ylabel('Disp (cm)',fontsize=14)
plt.xlim(xlim1)
plt.subplot(435)
plt.plot(de[0].times()-dt1,de[0].data*100)
plt.xlim(xlim1)
plt.subplot(436)
plt.plot(dz[0].times()-dt1,dz[0].data*100)
plt.xlim(xlim1)



plt.subplot(437)
plt.plot(n2[0].times()-dt1,n2[0].data*100,c='#FF8C00')
plt.xlim(xlim1)
plt.legend(['PLIG'],fontsize=14)
plt.ylabel('Accel (cm/s/s)',fontsize=14)
plt.subplot(438)
plt.plot(e2[0].times()-dt1,e2[0].data*100,c='#FF8C00')
plt.xlim(xlim1)
plt.subplot(439)
plt.plot(z2[0].times()-dt1,z2[0].data*100,c='#FF8C00')
plt.xlim(xlim1)

plt.subplot(4,3,10)
plt.plot(dn2[0].times()-dt1,dn2[0].data*100,c='#FF8C00')
plt.xlim(xlim1)
plt.ylabel('Disp (cm)',fontsize=14)
plt.xlabel('Seconds',fontsize=14)
plt.subplot(4,3,11)
plt.plot(de2[0].times()-dt1,de2[0].data*100,c='#FF8C00')
plt.xlim(xlim1)
plt.xlabel('Seconds',fontsize=14)
plt.subplot(4,3,12)
plt.plot(dz2[0].times()-dt1,dz2[0].data*100,c='#FF8C00')
plt.xlim(xlim1)
plt.xlabel('Seconds',fontsize=14)

plt.subplots_adjust(right=0.98,bottom=0.07)



plt.show()

