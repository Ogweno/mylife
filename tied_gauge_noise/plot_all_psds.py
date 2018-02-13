from mtspec import mtspec
from numpy import diff,where,c_,savetxt,genfromtxt,array,linspace,pi
from obspy import read
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

apgs=['M18B','M10B','G34B','FS20B','FS13B']
apg_depths=['731','685','2965','2389','2343']
darts=['46404','46407','46411','46419']
dart_depths=['2736','3268','4286','2773']
path_apgs='/Users/dmelgar/tidegauge_noise/APG/spectra/'
path_darts='/Users/dmelgar/tidegauge_noise/DART/spectra/'
months=['06','06','06','07']
xl=[12,1600]
yl=[1e-5,1e6]


xspec=linspace(12,120)
E0=1.5e-5
yspec=(E0/(8*pi))*xspec**2

plt.figure(figsize=(12,4))

plt.subplot(131)
plt.fill_between([12,120],yl[0],yl[1],color='b',alpha=0.2)
for ksta in range(len(apgs)):

    sta=apgs[ksta]

    spec=genfromtxt(path_apgs+sta+'_2017_05.spec')
    
    plt.loglog(spec[:,0],spec[:,1],label=sta+'('+apg_depths[ksta]+'m)')

plt.annotate(s='Tsunami band',xy=(13,10),fontsize=14,color='b')
plt.legend()
plt.grid()
plt.xlabel('Period (minutes)')
plt.ylabel('PSD')
plt.xlim(xl)
plt.ylim(yl)
plt.title('CI APGs')
plt.plot(xspec,yspec,'--',c='k',lw=2)
    
plt.subplot(132)
plt.fill_between([12,120],yl[0],yl[1],color='b',alpha=0.2)
for ksta in range(len(darts)):

    sta=darts[ksta]

    spec=genfromtxt(path_darts+sta+'_2017_'+months[ksta]+'.spec')
    
    plt.loglog(spec[:,0],spec[:,1],label=sta+'('+dart_depths[ksta]+'m)')        
                
plt.annotate(s='Tsunami band',xy=(13,10),fontsize=14,color='b')   
plt.legend()
plt.grid()
plt.xlabel('Period (minutes)')
plt.title('DARTs')
plt.xlim(xl)
plt.ylim(yl)
plt.plot(xspec,yspec,'--',c='k',lw=2)



plt.subplot(133)
plt.plot(1e-10,1e-10,c='#B22222',lw=1,label='DARTs')
plt.plot(1e-10,1e-10,c='#3CB371',lw=1,label='APGs')

plt.fill_between([12,120],yl[0],yl[1],color='b',alpha=0.2)
for ksta in range(len(darts)):

    sta=darts[ksta]
    spec=genfromtxt(path_darts+sta+'_2017_'+months[ksta]+'.spec')
    plt.loglog(spec[:,0],spec[:,1],c='#B22222',lw=0.3)        

for ksta in range(len(apgs)):
    sta=apgs[ksta]
    spec=genfromtxt(path_apgs+sta+'_2017_05.spec')
    plt.loglog(spec[:,0],spec[:,1],c='#3CB371',lw=0.3)             
    
plt.annotate(s='Tsunami band',xy=(13,10),fontsize=14,color='b')
plt.legend()
plt.grid()
plt.xlabel('Period (minutes)')
plt.title('All sites')
plt.xlim(xl)
plt.ylim(yl)
plt.plot(xspec,yspec,'--',c='k',lw=2)

plt.subplots_adjust(left=0.06,right=0.97,top=0.93,bottom=0.13)
plt.show() 