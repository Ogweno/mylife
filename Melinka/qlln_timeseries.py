from numpy import genfromtxt,c_,r_,ones,where,zeros,sin,pi,diff,arange,exp
from scipy.linalg import lstsq
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams

c1='#FF3333'
c2='#32CD32'
c3='#6495ED'

g=genfromtxt('/Users/dmelgar/Melinka2016/GPS/daily/QLLN.txt')
post_seismic=False
step_day=1089

# maximum day for linear fit
maxday=1089
epsilon=20

#get tiems
year=g[:,1]-2014
day=g[:,2]
t=day+year*365

#get data
eraw=g[:,3]
nraw=g[:,4]
uraw=g[:,5]
eerror=g[:,6]
nerror=g[:,7]
uerror=g[:,8]

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
def fitfunc(initial_guess,x,fix_step=None,ps_time=None,post_seismic=False):
    
    #Parse out arguemnts
    slope=initial_guess[0]
    intercept=initial_guess[1]
    step_time=initial_guess[2]
    step_before=initial_guess[3]
    step_after=initial_guess[4]
    yearly_amplitude=initial_guess[5]
    yearly_phase=initial_guess[6]
    semi_yearly_amplitude=initial_guess[7]
    semi_yearly_phase=initial_guess[8]
    ps_amplitude=initial_guess[9]
    decay=initial_guess[10]
    
    
    
    #initalize
    y = zeros(x.shape)
    
    ##Apply step
    if fix_step==None:
        y[x < step_time] = step_before
        y[step_time < x] = step_after
    else:
        y[x < fix_step] = step_before
        y[fix_step < x] = step_after        
    
    #Apply linear trend
    y=y+intercept+slope*x
    
    #Apply annual term
    y=y+yearly_amplitude*sin(2*pi*(1./365)*x+yearly_phase)
    
    #Apply semi-annual term
    y=y+semi_yearly_amplitude*sin(2*pi*(2./365)*x+semi_yearly_phase)
    
    #post-seismic term
    if post_seismic==True:
        i=where(x>fix_step)[0]
        y[i]=y[i]+ps_amplitude*exp((fix_step-x[i])/decay)
    
    return y


# Distance to the target function
errfunc = lambda p, x, y: fitfunc(p, x,fix_step=step_day,ps_time=step_day,post_seismic=post_seismic) - y 

#inital guess
slope=0.01
intercept=0
step_time=maxday
step_before=0
step_after=0
yearly_amplitude=5
yearly_phase=0
semi_yearly_amplitude=5
semi_yearly_phase=0
ps_time=maxday
ps_amplitude=1.0
decay=1.0
p0=[slope,intercept,step_time,step_before,step_after,yearly_amplitude,yearly_phase,
    semi_yearly_amplitude,semi_yearly_phase,ps_time,ps_amplitude,decay]

#get solutions
tfit=arange(t.min(),t.max()+1)

east_model, success = leastsq(errfunc, p0, args=(te, e))
yeast=fitfunc(east_model,tfit)

north_model, success = leastsq(errfunc, p0, args=(tn, n))
ynorth=fitfunc(north_model,tfit)

up_model, success = leastsq(errfunc, p0, args=(tu, u),epsfcn=10)
yup=fitfunc(up_model,tfit)





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
plt.plot(tfit,yeast,'k')
plt.title('East (mm)')
plt.xlim([t.min(),t.max()])
plt.ylim([-149,50])
plt.annotate('A',xy=(50,-125),fontsize=16)

ax=plt.subplot(332)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.scatter(tn,n,s=5,c=c2,lw=0)
plt.plot(tfit,ynorth,'k')
plt.title('North (mm)')
plt.xlim([t.min(),t.max()])
plt.ylim([-19,51])
plt.annotate('B',xy=(50,-10),fontsize=16)

ax=plt.subplot(333)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(4)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.scatter(tu,u,s=5,c=c3,lw=0)
plt.plot(tfit,yup,'k')
plt.title('Up (mm)')
plt.xlim([t.min(),t.max()])
plt.ylim([-58,58])
plt.annotate('C',xy=(50,-46),fontsize=16)



#get solutions
yeast2=fitfunc(east_model,te)
ynorth2=fitfunc(north_model,tn)
yup2=fitfunc(up_model,tu)



ax=plt.subplot(334)
plt.scatter(te,e-yeast2,s=5,c=c1,lw=0)
#plt.plot(tfit,yeast,'k')
plt.xlim([t.min(),t.max()])
ymajorLocator = MultipleLocator(10)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.annotate('D',xy=(50,-12),fontsize=16)

ax=plt.subplot(335)
plt.scatter(tn,n-ynorth2,s=5,c=c2,lw=0)
#plt.plot(tfit,ynorth,'k')
plt.xlim([t.min(),t.max()])
ymajorLocator = MultipleLocator(4)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.annotate('E',xy=(50,-4.5),fontsize=16)

ax=plt.subplot(336)
plt.scatter(tu,u-yup2,s=5,c=c3,lw=0)
#plt.plot(tfit,yup,'k')
plt.xlim([t.min(),t.max()])
ymajorLocator = MultipleLocator(10)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.annotate('F',xy=(50,-14),fontsize=16)



xmajorLocator = MultipleLocator(50)
xminorLocator = MultipleLocator(5)

ax=plt.subplot(337)
plt.scatter(te,e,s=8,c=c1,lw=0)
plt.plot(tfit,yeast,'k')
plt.xlim([950,t.max()])
ymajorLocator = MultipleLocator(50)
yminorLocator = MultipleLocator(10)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.annotate('G',xy=(959,-127),fontsize=16)

ax=plt.subplot(338)
plt.scatter(tn,n,s=8,c=c2,lw=0)
plt.plot(tfit,ynorth,'k')
plt.xlabel('Days since Jan 1st 2014',fontsize=14)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.annotate('H',xy=(959,-10.67),fontsize=16)

plt.xlim([950,t.max()])

ax=plt.subplot(339)
plt.scatter(tu,u,s=8,c=c3,lw=0)
plt.plot(tfit,yup,'k')
plt.xlim([950,t.max()])
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(4)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
plt.annotate('I',xy=(959,-46),fontsize=16)

plt.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.1)
plt.show()