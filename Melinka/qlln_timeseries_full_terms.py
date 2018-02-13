from numpy import genfromtxt,c_,r_,ones,where,zeros,sin,pi,diff,arange,exp
from scipy.linalg import lstsq
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams
from obspy.core import UTCDateTime

c1='#FF3333'
c2='#32CD32'
c3='#6495ED'

#Read new data
date=genfromtxt('/Users/dmelgar/Melinka2016/GPS/daily/qlln.crd',usecols=1,dtype='S')
g=genfromtxt('/Users/dmelgar/Melinka2016/GPS/daily/qlln.crd')
eraw=g[:,4]*1000
nraw=g[:,3]*1000
uraw=g[:,5]*1000
eraw=eraw-eraw[0]
nraw=nraw-nraw[0]
uraw=uraw-uraw[0]
eerror=g[:,7]
nerror=g[:,6]
uerror=g[:,8]
t0=UTCDateTime(date[0])
t=zeros(len(date))
for k in range(len(date)):
    t[k]=(UTCDateTime(date[k])-t0)/86400
#New run time params
post_seismic=True
step_day=724
ps_start=725
epsilon=20
ps_amplitude=10
decay=20





#remove outliers and mid step point
def despike(t,y,error,thresh,maxt):
    i=where(t<maxt)[0]
    keep=where(t>=maxt)[0]
    j=where(abs(diff(y[i]))<thresh)[0]
    final=r_[i[j],keep]
    tout=t[final]
    yout=y[final]
    errout=error[final]
    return tout,yout,errout

te,e,ee=despike(t,eraw,eerror,12.0,t.max())  
tn,n,ne=despike(t,nraw,nerror,5.0,t.max())  
tu,u,ue=despike(t,uraw,uerror,22,t.max()) 

#remove midstep day
i=where(te!=1090)[0]
te=te[i]
e=e[i]
eerror=ee[i]
i=where(tn!=1090)[0]
tn=tn[i]
n=n[i]
nerror=ne[i]
i=where(tu!=1090)[0]
tu=tu[i]
u=u[i]
uerror=ue[i]

#Functional form of time series
#Functional form of time series
def fitfunc(initial_guess,x,step_day,ps_start):
    
    #Parse out arguemnts
    slope=initial_guess[0]
    intercept=initial_guess[1]
    step_before=initial_guess[2]
    step_after=initial_guess[3]
    yearly_amplitude=initial_guess[4]
    yearly_phase=initial_guess[5]
    semi_yearly_amplitude=initial_guess[6]
    semi_yearly_phase=initial_guess[7]
    ps_amplitude=initial_guess[8]
    decay=initial_guess[9]
    
    #initalize
    y = zeros(x.shape)
    
    ##Apply step
    y[x <= step_day] = step_before
    y[step_day < x] = step_after        
    
    #Apply linear trend
    y=y+intercept+slope*x
    
    #Apply annual term
    y=y+yearly_amplitude*sin(2*pi*(1./365)*x+yearly_phase)
    
    #Apply semi-annual term
    y=y+semi_yearly_amplitude*sin(2*pi*(2./365)*x+semi_yearly_phase)
    
    #post seismic
    i=where(x >= ps_start)[0]
    y[i]=y[i]+ps_amplitude*exp((x[i]-ps_start)/decay)

    
    return y
    
def clean(data,model,x,step_day,ps_start):
    
    #Parse out arguemnts
    slope=model[0]
    intercept=model[1]
    step_before=model[2]
    step_after=model[3]
    yearly_amplitude=model[4]
    yearly_phase=model[5]
    semi_yearly_amplitude=model[6]
    semi_yearly_phase=model[7]
    ps_amplitude=model[8]
    decay=model[9]     
    
    #Apply linear trend
    data=data-slope*x
    
    #Apply annual term
    data=data-yearly_amplitude*sin(2*pi*(1./365)*x+yearly_phase)
    
    #Apply semi-annual term
    data=data-semi_yearly_amplitude*sin(2*pi*(2./365)*x+semi_yearly_phase)

    return data

# Distance to the target function
errfunc = lambda p, x, y: fitfunc(p, x,step_day,ps_start) - y 

#inital guess
slope=0.01
intercept=0
step_before=0
step_after=-70
yearly_amplitude=5
yearly_phase=0
semi_yearly_amplitude=5
semi_yearly_phase=0
eps_amplitude=-10
nps_amplitude=-5
ups_amplitude=-5
ps_decay=10

p0=[slope,intercept,step_before,step_after,yearly_amplitude,yearly_phase,
    semi_yearly_amplitude,semi_yearly_phase,eps_amplitude,ps_decay]
east_model, success = leastsq(errfunc, p0, args=(te, e))
yeast=fitfunc(east_model,te,step_day,ps_start)

p0=[slope,intercept,step_before,step_after,yearly_amplitude,yearly_phase,
    semi_yearly_amplitude,semi_yearly_phase,nps_amplitude,ps_decay]
north_model, success = leastsq(errfunc, p0, args=(tn, n))
ynorth=fitfunc(north_model,tn,step_day,ps_start)

p0=[slope,intercept,step_before,step_after,yearly_amplitude,yearly_phase,
    semi_yearly_amplitude,semi_yearly_phase,ups_amplitude,ps_decay]
up_model, success = leastsq(errfunc, p0, args=(tu, u))
yup=fitfunc(up_model,tu,step_day,ps_start)


#remove seasonal,interseismic, and coseismic terms
ye_clean=clean(e,east_model,te,step_day,ps_start)
yn_clean=clean(n,north_model,tn,step_day,ps_start)
yu_clean=clean(u,up_model,tu,step_day,ps_start)


# Remove everything
ye_full_clean=e-yeast
yn_full_clean=n-ynorth
yu_full_clean=u-yup




plt.figure(figsize=(14,7))

rcParams['xtick.major.size'] = 6.0
rcParams['xtick.major.width'] = 0.5
rcParams['xtick.minor.size'] = 4
rcParams['xtick.minor.width'] = 0.5
rcParams['ytick.major.size'] = 6.0
rcParams['ytick.major.width'] = 0.5
rcParams['ytick.minor.size'] = 4
rcParams['ytick.minor.width'] = 0.5
xmajorLocator = MultipleLocator(200)
xminorLocator = MultipleLocator(20)

ax=plt.subplot(331)
ymajorLocator = MultipleLocator(50)
yminorLocator = MultipleLocator(10)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(te,e,s=5,c=c1,lw=0)
ax.plot(te,yeast,'k')
ax.set_title('East (mm)')
ax.set_xlim([t.min(),t.max()])
ax.set_ylim([-170,40])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('A',xy=(50,-145),fontsize=16)

ax=plt.subplot(332)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(tn,n,s=5,c=c2,lw=0)
ax.plot(tn,ynorth,'k')
ax.set_title('North (mm)')
ax.set_xlim([t.min(),t.max()])
ax.set_ylim([-36,32])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('B',xy=(50,-24),fontsize=16)

ax=plt.subplot(333)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(4)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(tu,u,s=5,c=c3,lw=0)
ax.plot(tu,yup,'k')
ax.set_title('Up (mm)')
ax.set_xlim([t.min(),t.max()])
ax.set_ylim([-58,58])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('C',xy=(50,-40),fontsize=16)







ax=plt.subplot(334)
ymajorLocator = MultipleLocator(50)
yminorLocator = MultipleLocator(10)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(te,ye_clean,s=5,c=c1,lw=0)
ax.set_xlim([t.min(),t.max()])
ax.set_ylim([-200,30])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('D',xy=(50,-165),fontsize=16)

ax=plt.subplot(335)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(tn,yn_clean,s=5,c=c2,lw=0)
ax.set_xlim([t.min(),t.max()])
ax.set_ylim([-70,10])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('E',xy=(50,-54),fontsize=16)

ax=plt.subplot(336)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(4)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(tu,yu_clean,s=5,c=c3,lw=0)
ax.set_xlim([t.min(),t.max()])
ax.set_ylim([-70,25])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('F',xy=(50,-54),fontsize=16)





xmajorLocator = MultipleLocator(20)
xminorLocator = MultipleLocator(4)

ax=plt.subplot(337)
ymajorLocator = MultipleLocator(5)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(te,ye_full_clean,s=8,c=c1,lw=0)
ax.plot([724,724],[-20,20],'--',c='#606060',lw=2)
ax.plot([0,852],[0,0],c='#606060',lw=1)
ax.set_xlim([724-90,724+90])
ax.set_ylim([-11,11])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('G',xy=(645,-7),fontsize=16)

ax=plt.subplot(338)
ymajorLocator = MultipleLocator(5)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(tn,yn_full_clean,s=8,c=c2,lw=0)
ax.plot([724,724],[-20,20],'--',c='#606060',lw=2)
ax.plot([0,852],[0,0],c='#606060',lw=1)
ax.set_xlim([724-90,724+90])
ax.set_ylim([-11,11])
ax.set_xlabel('Days since Jan 1st 2015',fontsize=14)
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('H',xy=(645,-7),fontsize=16)


ax=plt.subplot(339)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.scatter(tu,yu_full_clean,s=8,c=c3,lw=0)
ax.plot([724,724],[-60,60],'--',c='#606060',lw=2)
ax.plot([0,852],[0,0],c='#606060',lw=1)
ax.set_xlim([724-90,724+90])
ax.set_ylim([-41,41])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('I',xy=(645,-27),fontsize=16)

plt.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.1)



plt.show()
