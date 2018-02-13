from matplotlib import pyplot as plt
import gmpe_tools
import bbptools
from numpy import logspace,zeros,ones,exp,log,log10,array

path='/Users/dmelgar/code/BBP/bbp/bbp_data/finished/eew_test/'
sims=['45','50','55','60','65','70']
xl=[1,50]

### PGA
plt.figure(figsize=(13,8))


#######      M4.5     ##########
k=0
sim_path=path+sims[k]

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgalims=log10(pga)
pga_range=pgalims.max()-pgalims.min()
yl=[0.001,0.6]

#get predicted data
Rjb_pred=logspace(-1,2)
M=M*ones(len(Rjb_pred))
vs30=720*ones(len(Rjb_pred))
U=zeros(len(Rjb_pred))
NS=zeros(len(Rjb_pred))
RS=zeros(len(Rjb_pred))
SS=ones(len(Rjb_pred))
pred_pga,std_pga=gmpe_tools.bssa14(M,Rjb_pred,vs30,SS=SS,NS=NS,RS=RS,U=U,intensity_measure='PGA')

plt.subplot(231)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(Rjb,pga,lw=0.5,s=40,c='#DC143C')
plt.plot(Rjb_pred,pred_pga,lw=2,c='#DC143C')
plt.plot(Rjb_pred,exp(log(pred_pga)+1.65*std_pga),'--',lw=2,c='#DC143C')
plt.plot(Rjb_pred,exp(log(pred_pga)-1.65*std_pga),'--',lw=2,c='#DC143C')
plt.xlim(xl)
plt.ylim(yl)
plt.ylabel('PGA g')
plt.legend(['M4.5'],loc=3,frameon=False)

########       END SUBPLOT     ########


#######      M5.0     ##########
k=1
sim_path=path+sims[k]

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgalims=log10(pga)
pga_range=pgalims.max()-pgalims.min()
yl=[0.002,1.1]

#get predicted data
Rjb_pred=logspace(-1,2)
M=M*ones(len(Rjb_pred))
vs30=720*ones(len(Rjb_pred))
U=zeros(len(Rjb_pred))
NS=zeros(len(Rjb_pred))
RS=zeros(len(Rjb_pred))
SS=ones(len(Rjb_pred))
pred_pga,std_pga=gmpe_tools.bssa14(M,Rjb_pred,vs30,SS=SS,NS=NS,RS=RS,U=U,intensity_measure='PGA')

plt.subplot(232)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(Rjb,pga,lw=0.5,s=40,c='#32CD32')
plt.plot(Rjb_pred,pred_pga,lw=2,c='#32CD32')
plt.plot(Rjb_pred,exp(log(pred_pga)+1.65*std_pga),'--',lw=2,c='#32CD32')
plt.plot(Rjb_pred,exp(log(pred_pga)-1.65*std_pga),'--',lw=2,c='#32CD32')
plt.xlim(xl)
plt.ylim(yl)
plt.legend(['M5.0'],loc=3,frameon=False)

########       END SUBPLOT     ########


#######      M5.5     ##########
k=2
sim_path=path+sims[k]

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgalims=log10(pga)
pga_range=pgalims.max()-pgalims.min()
yl=[0.008,1.7]

#get predicted data
Rjb_pred=logspace(-1,2)
M=M*ones(len(Rjb_pred))
vs30=720*ones(len(Rjb_pred))
U=zeros(len(Rjb_pred))
NS=zeros(len(Rjb_pred))
RS=zeros(len(Rjb_pred))
SS=ones(len(Rjb_pred))
pred_pga,std_pga=gmpe_tools.bssa14(M,Rjb_pred,vs30,SS=SS,NS=NS,RS=RS,U=U,intensity_measure='PGA')

plt.subplot(233)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(Rjb,pga,lw=0.5,s=40,c='#0000CD')
plt.plot(Rjb_pred,pred_pga,lw=2,c='#0000CD')
plt.plot(Rjb_pred,exp(log(pred_pga)+1.65*std_pga),'--',lw=2,c='#0000CD')
plt.plot(Rjb_pred,exp(log(pred_pga)-1.65*std_pga),'--',lw=2,c='#0000CD')
plt.xlim(xl)
plt.ylim(yl)
plt.legend(['M5.5'],loc=3,frameon=False)

########       END SUBPLOT     ########


#######      M6.0     ##########
k=3
sim_path=path+sims[k]

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgalims=log10(pga)
pga_range=pgalims.max()-pgalims.min()
yl=[0.01,1.2]

#get predicted data
Rjb_pred=logspace(-1,2)
M=M*ones(len(Rjb_pred))
vs30=720*ones(len(Rjb_pred))
U=zeros(len(Rjb_pred))
NS=zeros(len(Rjb_pred))
RS=zeros(len(Rjb_pred))
SS=ones(len(Rjb_pred))
pred_pga,std_pga=gmpe_tools.bssa14(M,Rjb_pred,vs30,SS=SS,NS=NS,RS=RS,U=U,intensity_measure='PGA')

plt.subplot(234)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(Rjb,pga,lw=0.5,s=40,c='#DAA520')
plt.plot(Rjb_pred,pred_pga,lw=2,c='#DAA520')
plt.plot(Rjb_pred,exp(log(pred_pga)+1.65*std_pga),'--',lw=2,c='#DAA520')
plt.plot(Rjb_pred,exp(log(pred_pga)-1.65*std_pga),'--',lw=2,c='#DAA520')
plt.xlim(xl)
plt.ylim(yl)
plt.ylabel('PGA g')
plt.xlabel('Rjb (km)')
plt.legend(['M6.0'],loc=3,frameon=False)

########       END SUBPLOT     ########



#######      M6.5     ##########
k=4
sim_path=path+sims[k]

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgalims=log10(pga)
pga_range=pgalims.max()-pgalims.min()
yl=[0.01,1.2]

#get predicted data
Rjb_pred=logspace(-1,2)
M=M*ones(len(Rjb_pred))
vs30=720*ones(len(Rjb_pred))
U=zeros(len(Rjb_pred))
NS=zeros(len(Rjb_pred))
RS=zeros(len(Rjb_pred))
SS=ones(len(Rjb_pred))
pred_pga,std_pga=gmpe_tools.bssa14(M,Rjb_pred,vs30,SS=SS,NS=NS,RS=RS,U=U,intensity_measure='PGA')

plt.subplot(235)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(Rjb,pga,lw=0.5,s=40,c='#9932CC')
plt.plot(Rjb_pred,pred_pga,lw=2,c='#9932CC')
plt.plot(Rjb_pred,exp(log(pred_pga)+1.65*std_pga),'--',lw=2,c='#9932CC')
plt.plot(Rjb_pred,exp(log(pred_pga)-1.65*std_pga),'--',lw=2,c='#9932CC')
plt.xlim(xl)
plt.ylim(yl)
plt.xlabel('Rjb (km)')
plt.legend(['M6.5'],loc=3,frameon=False)

########       END SUBPLOT     ########



#######      M7.0    ##########
k=5
sim_path=path+sims[k]

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgalims=log10(pga)
pga_range=pgalims.max()-pgalims.min()
yl=[0.02,1.4]
#get predicted data
Rjb_pred=logspace(-1,2)
M=M*ones(len(Rjb_pred))
vs30=720*ones(len(Rjb_pred))
U=zeros(len(Rjb_pred))
NS=zeros(len(Rjb_pred))
RS=zeros(len(Rjb_pred))
SS=ones(len(Rjb_pred))
pred_pga,std_pga=gmpe_tools.bssa14(M,Rjb_pred,vs30,SS=SS,NS=NS,RS=RS,U=U,intensity_measure='PGA')

plt.subplot(236)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(Rjb,pga,lw=0.5,s=40,c='#202020')
plt.plot(Rjb_pred,pred_pga,lw=2,c='#202020')
plt.plot(Rjb_pred,exp(log(pred_pga)+1.65*std_pga),'--',lw=2,c='#202020')
plt.plot(Rjb_pred,exp(log(pred_pga)-1.65*std_pga),'--',lw=2,c='#202020')
plt.xlim(xl)
plt.ylim(yl)
plt.xlabel('Rjb (km)')
plt.legend(['M7.0'],loc=3,frameon=False)

########       END SUBPLOT     ########

plt.subplots_adjust(left=0.09,right=0.96,bottom=0.09,top=0.96,wspace=0.12,hspace=0.09)

plt.show()