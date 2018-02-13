from numpy import genfromtxt,where,arange,argmax
from matplotlib import pyplot as plt
from scipy.signal import correlate
from scipy.interpolate import interp1d

static=genfromtxt(u'/Users/dmelgar/Tsunamis/maule_gps_tg_static/_output/fort.gauge')
kinematic=genfromtxt(u'/Users/dmelgar/Tsunamis/maule_gps_tg_kinematic/_output/fort.gauge')
fout="/Users/dmelgar/Maule2010/plots/tsunami_profiles/dip4.pdf"
profile=range(70,80)

#Init
fig, axarr = plt.subplots(10, 1,figsize=(6, 11))  
fig.set_figwidth(5)
fig.set_figheight(10)
k=0
xl=[0,90]
annotx=70
for sta in profile:
    i_st=where(static[:,0]==sta)[0]
    i_ki=where(kinematic[:,0]==sta)[0]
    t_st=static[i_st,2]
    eta_st=static[i_st,6]
    t_ki=kinematic[i_ki,2]
    eta_ki=kinematic[i_ki,6]
    #Resample to regular interval
    ti=arange(0,7199,1)
    f_st=interp1d(t_st,eta_st)
    f_ki=interp1d(t_ki,eta_ki)
    eta_st_i=f_st(ti)
    eta_ki_i=f_ki(ti)
    #get correlation lag
    c=correlate(eta_st_i,eta_ki_i)
    lag=abs(argmax(c)-len(c)/2)
    #Make plots
    ax=axarr[k]
    ax.plot(ti*(1./60),eta_st_i,'k',label='Instantaneous')
    ax.plot(ti*(1./60),eta_ki_i,'r',label='Kinematic')
    ax.set_ylabel(r'$\eta$ (m)')
    if k==9:
        ax.set_xlabel('Minutes after OT')
    if k==0:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.8),
          ncol=2, fancybox=True, shadow=True)
    if k!=9:
        ax.xaxis.set_ticklabels([])
    k+=1
    ax.set_xlim(xl)
    ax.set_ylim([1.1*eta_st.min(),1.1*eta_st.max()])
    ax.annotate(r'$\tau$ = '+str(lag)+'s',xy=(annotx,0.6*eta_st.max()))
    ax.locator_params(axis='y',nbins=4)
plt.subplots_adjust(hspace=0.2)
plt.savefig(fout)
plt.show()