from matplotlib import pyplot as plt
import gmpe_tools
import bbptools
from numpy import logspace,zeros,ones,exp,log,log10,array,r_,arange

magnitude='m7.0'
M=7.0
path='/Users/dmelgar/code/BBP/bbp/bbp_data/finished/large_eew_test/'
sims=['0.05','0.25','0.5','0.75','0.95']
xl=[1,50]
yl=[0.001,0.6]
yl=[0.019,1.5]


### PGA



PGA=array([])
PGV=array([])
PGD=array([])
Rjb=array([])

#######      M4.5     ##########
k=0
sim_path=path+sims[k]

for k in range(len(sims)):
    #get simulation data
    sim_path=path+magnitude+'_frac'+sims[k]
    pga,pgv,pgd,rjb,foo=bbptools.get_simulation_intensity(sim_path)
    PGA=r_[PGA,pga]
    PGV=r_[PGV,pgv]
    PGD=r_[PGD,pgd]
    Rjb=r_[Rjb,rjb]
    
pgalims=log10(PGA)
pga_range=pgalims.max()-pgalims.min()


plt.figure(figsize=(5,5))

#get predicted data
Rjb_pred=logspace(-1,2)
M=M*ones(len(Rjb_pred))
vs30=720*ones(len(Rjb_pred))
U=zeros(len(Rjb_pred))
NS=zeros(len(Rjb_pred))
RS=zeros(len(Rjb_pred))
SS=ones(len(Rjb_pred))
pred_pga,std_pga=gmpe_tools.bssa14(M,Rjb_pred,vs30,SS=SS,NS=NS,RS=RS,U=U,intensity_measure='PGA')

ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')

k=0
i=arange(46*k,46*(k+1))
plt.scatter(Rjb[i],PGA[i],lw=0.5,s=40,c='#DC143C')
k=1
i=arange(46*k,46*(k+1))
plt.scatter(Rjb[i],PGA[i],lw=0.5,s=40,c='#32CD32')
k=2
i=arange(46*k,46*(k+1))
plt.scatter(Rjb[i],PGA[i],lw=0.5,s=40,c='#0000CD')
k=3
i=arange(46*k,46*(k+1))
plt.scatter(Rjb[i],PGA[i],lw=0.5,s=40,c='#DAA520')
k=4
i=arange(46*k,46*(k+1))
plt.scatter(Rjb[i],PGA[i],lw=0.5,s=40,c='#9932CC')

plt.legend(['5%','25%','50%','75%','95%'],loc=3,frameon=False)

plt.plot(Rjb_pred,pred_pga,lw=2,c='#404040')
plt.plot(Rjb_pred,exp(log(pred_pga)+1.65*std_pga),'--',lw=2,c='#404040')
plt.plot(Rjb_pred,exp(log(pred_pga)-1.65*std_pga),'--',lw=2,c='#404040')
plt.xlim(xl)
plt.ylim(yl)
plt.ylabel('PGA g')
plt.title(magnitude.upper())



plt.subplots_adjust(left=0.15,right=0.97,top=0.94,bottom=0.11)

plt.show()