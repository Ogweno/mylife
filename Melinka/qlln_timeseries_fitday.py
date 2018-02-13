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

#Read old data
#g=genfromtxt('/Users/dmelgar/Melinka2016/GPS/daily/QLLN.txt')
##get times
#year=g[:,1]-2014
#day=g[:,2]
#t=day+year*365
#eraw=g[:,3]
#nraw=g[:,4]
#uraw=g[:,5]
#eerror=g[:,6]
#nerror=g[:,7]
#uerror=g[:,8]
##Old run time params
#post_seismic=False
#step_day=1089
#fit_day=1087
#maxday=1089
#epsilon=20


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
post_seismic=False
step_day=724
ps_start=725
fit_day=1000 #fit to this day
maxday=724
epsilon=20






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

te,e,ee=despike(t,eraw,eerror,12.0,maxday)  
tn,n,ne=despike(t,nraw,nerror,5.0,maxday)  
tu,u,ue=despike(t,uraw,uerror,22,maxday) 

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
def fitfunc1(initial_guess,x,post_seismisc=False,ps_start=0):
    
    #Parse out arguemnts
    slope=initial_guess[0]
    intercept=initial_guess[1]
    yearly_amplitude=initial_guess[2]
    yearly_phase=initial_guess[3]
    semi_yearly_amplitude=initial_guess[4]
    semi_yearly_phase=initial_guess[5]
    ps_amplitude=initial_guess[6]
    decay=initial_guess[7]
    
    
    #initalize
    y = zeros(x.shape)
       
    
    #Apply linear trend
    y=y+intercept+slope*x
    
    #Apply annual term
    y=y+yearly_amplitude*sin(2*pi*(1./365)*x+yearly_phase)
    
    #Apply semi-annual term
    y=y+semi_yearly_amplitude*sin(2*pi*(2./365)*x+semi_yearly_phase)
    
    #post-seismic term
    if post_seismic==True:
        i=where(x>ps_start)[0]
        y[i]=y[i]+ps_amplitude*exp((ps_start-x[i])/decay)
    
    return y
    
    
def fitfunc2(initial_guess,x,fix_step=None):
    
    #Parse out arguemnts
    step_time=initial_guess[0]
    step_before=initial_guess[1]
    step_after=initial_guess[2]
    
    
    
    #initalize
    y = zeros(x.shape)
    
    ##Apply step
    if fix_step==None:
        y[x < step_time] = step_before
        y[step_time < x] = step_after
    else:
        y[x < fix_step] = step_before
        y[fix_step < x] = step_after        
    
    return y
    


# Distance to the target function
errfunc = lambda p, x, y: fitfunc1(p, x,post_seismisc=post_seismic,ps_start=ps_start) - y 
errfunc2 = lambda p, x, y: fitfunc2(p, x) - y 

#inital guess
slope=0.01
intercept=0
yearly_amplitude=5
yearly_phase=0
semi_yearly_amplitude=5
semi_yearly_phase=0
ps_amplitude=10
ps_decay=1
p0=[slope,intercept,yearly_amplitude,yearly_phase,
    semi_yearly_amplitude,semi_yearly_phase,ps_amplitude,ps_decay]


#Only fit to beofre day of earthquake
i=where(te<=fit_day)
tefit=te[i]
efit=e[i]
i=where(tn<=fit_day)
tnfit=tn[i]
nfit=n[i]
i=where(tu<=fit_day)
tufit=tu[i]
ufit=u[i]

#get solutions
tfit=arange(t.min(),t.max()+1)



east_model, success = leastsq(errfunc, p0, args=(tefit, efit))
yeast=fitfunc1(east_model,tefit)

north_model, success = leastsq(errfunc, p0, args=(tnfit, nfit))
ynorth=fitfunc1(north_model,tnfit)

up_model, success = leastsq(errfunc, p0, args=(tufit, ufit),epsfcn=10)
yup=fitfunc1(up_model,tufit)





# MODEL 2

step_time=maxday
step_before=0
step_after=0

p0=[step_time,step_before,step_after]

#get solutions, first remove seasonal terms
ye_season=fitfunc1(east_model,te)
yn_season=fitfunc1(north_model,tn)
yu_season=fitfunc1(up_model,tu)

ecorrect=e-ye_season
ncorrect=n-yn_season
ucorrect=u-yu_season

east_model2, success = leastsq(errfunc2, p0, args=(te, ecorrect))
north_model2, success = leastsq(errfunc2, p0, args=(tn, ncorrect))
up_model2, success = leastsq(errfunc2, p0, args=(tu, ucorrect),epsfcn=10)







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
plt.scatter(te,e,s=5,c=c1,lw=0)
plt.plot(tefit,yeast,'k')
plt.title('East (mm)')
plt.xlim([t.min(),t.max()])
plt.ylim([-149,50])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
plt.annotate('A',xy=(50,-125),fontsize=16)

ax=plt.subplot(332)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.scatter(tn,n,s=5,c=c2,lw=0)
plt.plot(tnfit,ynorth,'k')
plt.title('North (mm)')
plt.xlim([t.min(),t.max()])
plt.ylim([-19,51])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
plt.annotate('B',xy=(50,-10),fontsize=16)

ax=plt.subplot(333)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(4)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.scatter(tu,u,s=5,c=c3,lw=0)
plt.plot(tufit,yup,'k')
plt.title('Up (mm)')
plt.xlim([t.min(),t.max()])
plt.ylim([-58,58])
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
plt.annotate('C',xy=(50,-46),fontsize=16)



#subtract previous model
yeast2=fitfunc2(east_model2,tfit)
ynorth2=fitfunc2(north_model2,tfit)
yup2=fitfunc2(up_model2,tfit)



ax=plt.subplot(334)
plt.scatter(te,ecorrect,s=5,c=c1,lw=0)
plt.plot(tfit,yeast2,'k')
plt.xlim([t.min(),t.max()])
ymajorLocator = MultipleLocator(50)
yminorLocator = MultipleLocator(10)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
plt.annotate('D',xy=(50,-164),fontsize=16)

ax=plt.subplot(335)
plt.scatter(tn,ncorrect,s=5,c=c2,lw=0)
plt.plot(tfit,ynorth2,'k')
plt.xlim([t.min(),t.max()])
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
plt.annotate('E',xy=(50,-59),fontsize=16)

ax=plt.subplot(336)
plt.scatter(tu,ucorrect,s=5,c=c3,lw=0)
plt.plot(tfit,yup2,'k')
plt.xlim([t.min(),t.max()])
ymajorLocator = MultipleLocator(40)
yminorLocator = MultipleLocator(5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
plt.annotate('F',xy=(50,-79),fontsize=16)




yeast3=fitfunc2(east_model2,te)
ynorth3=fitfunc2(north_model2,tn)
yup3=fitfunc2(up_model2,tu)


xmajorLocator = MultipleLocator(20)
xminorLocator = MultipleLocator(4)

ax=plt.subplot(337)
ax.scatter(te,ecorrect-yeast3,s=8,c=c1,lw=0)
#plt.plot(tfit,yeast2,'k')
ax.plot([1089,1089],[-20,20],'--',c='#606060',lw=2)
ax.plot([0,2000],[0,0],c='#606060',lw=1)
ax.set_xlim([650,t.max()])
ax.set_ylim([-11,11])
ymajorLocator = MultipleLocator(5)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('G',xy=(1005,-7),fontsize=16)

ax=plt.subplot(338)
ax.scatter(tn,ncorrect-ynorth3,s=8,c=c2,lw=0)
ax.plot([1089,1089],[-20,20],'--',c='#606060',lw=2)
ax.plot([0,2000],[0,0],c='#606060',lw=1)
ax.set_xlim([650,t.max()])
ax.set_ylim([-11,11])
#ax.plot(tfit,ynorth2,'k')
ax.set_xlabel('Days since Jan 1st 2014',fontsize=14)
ymajorLocator = MultipleLocator(5)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('H',xy=(1005,-7),fontsize=16)


ax=plt.subplot(339)
ax.scatter(tu,ucorrect-yup3,s=8,c=c3,lw=0)
ax.plot([1089,1089],[-60,60],'--',c='#606060',lw=2)
ax.plot([0,2000],[0,0],c='#606060',lw=1)
ax.set_xlim([650,t.max()])
ax.set_ylim([-41,41])
#ax.plot(tfit,yup2,'k')
ax.set_xlim([650,t.max()])
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.get_xaxis().set_tick_params(which='both',direction='out')
ax.get_yaxis().set_tick_params(which='both',direction='out')
ax.annotate('I',xy=(1005,-27),fontsize=16)

plt.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.1)
plt.show()