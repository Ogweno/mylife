from numpy import genfromtxt,zeros,log10,arange,deg2rad,where
import matplotlib.path as mplPath
from scipy.stats import binned_statistic
from matplotlib import pyplot as plt
from matplotlib import rcParams

bin_edges=arange(-47,-35.9,0.2)
w=0.2
#30 arcsecond area
s=6371e3*deg2rad((30./3600))
area=s**2
mu=70e9

rcParams['xtick.major.size'] = 6.0
rcParams['xtick.major.width'] = 0.5
rcParams['xtick.minor.size'] = 4
rcParams['xtick.minor.width'] = 0.5

#read in gCMT
gcmt=genfromtxt('/Users/dmelgar/code/GMT/Melinka/gcmt_thrust.psmeca',usecols=[0,1,3,4,5,6,7,8,9])
lon_gcmt=gcmt[:,0]
lon_gcmt=gcmt[:,1]

#selection polygon
poly=genfromtxt(u'/Users/dmelgar/code/GMT/Melinka/trench_poly.txt')

#select points in polygon
bbPath = mplPath.Path(poly)
i=bbPath.contains_points(gcmt[:,0:2])
events=gcmt[i,:]

moment=zeros(len(events))
Mw=zeros(len(events))
for k in range(len(events)):
    moment[k]=(1./2**0.5)*(10**events[k,8]*events[k,2]**2+10**events[k,8]*events[k,3]**2+10**events[k,8]*events[k,4]**2+10**events[k,8]*events[k,5]**2+10**events[k,8]*events[k,6]**2+10**events[k,7]*events[k,2]**2)/1e7
    Mw[k]=(2./3)*(log10(moment[k])-9.1)
    
#bin by latitude
bin_moment,a,b=binned_statistic(events[:,1], moment, statistic='sum', bins=bin_edges, range=None)
i=where(bin_moment<1e12)[0]
bin_moment[i]=1e12


#load slip ivnersion and get binend moment
fault=genfromtxt('/Users/dmelgar/Slip_inv/Melinka_usgs/output/inverse_models/models/gps_sm_insar_w3.01.0011.inv.total')
fault_moment=fault[:,10]*fault[:,11]*fault[:,12]*(fault[:,8]**2+fault[:,9]**2)**0.5
bin_fault_moment,a,b=binned_statistic(fault[:,2], fault_moment, statistic='sum', bins=bin_edges, range=None)
i=where(bin_fault_moment<1e12)[0]
bin_fault_moment[i]=1e12
    
#load valdivia
val=genfromtxt(u'/Users/dmelgar/code/GMT/Melinka/valdivia.xyz')
val_moment=area*val[:,2]*mu
val_Mw=(2./3)*(log10(val_moment.sum())-9.1)
print 'Valdivia magnitude is '+str(val_Mw)
bin_val_moment,a,b=binned_statistic(val[:,1], val_moment, statistic='sum', bins=bin_edges, range=None)
i=where(bin_val_moment<1e12)[0]
bin_val_moment[i]=1e12

#load maule
maul=genfromtxt('/Users/dmelgar/Slip_inv/maule_gps_tg/output/inverse_models/models/gps_tg_10win3.5_tgdiv2_nocorr.0003.inv.total')
maul_moment=maul[:,10]*maul[:,11]*maul[:,12]*(maul[:,8]**2+maul[:,9]**2)**0.5
bin_maul_moment,a,b=binned_statistic(maul[:,2], maul_moment, statistic='sum', bins=bin_edges, range=None)
i=where(bin_maul_moment<1e12)[0]
bin_maul_moment[i]=1e12

#plot
fig=plt.figure(figsize=(3.0,14))
ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')

#pplot stuff for legend
plt.barh(bin_edges[0],1e12,color='#228B22',height=w)
plt.barh(bin_edges[0],1e12,color='#FFC125',height=w)
plt.barh(bin_edges[0],1e12,color='#EE2C2C',height=w)
plt.barh(bin_edges[0],1e12,color='#3D59AB',height=w)
plt.legend(['1960 M9.5 Valdivia','2010 M8.8 Maule','2016 M7.6 Melinka','Events 1976-2016'],fontsize=13,bbox_to_anchor=(1.04, 1.15),frameon=False)

#now go one bin at a time and plot them in ascending order
for k in range(len(bin_edges)-1):
    if bin_val_moment[k]>bin_maul_moment[k]:
        plt.barh(bin_edges[k],bin_val_moment[k],color='#228B22',height=w)
        plt.barh(bin_edges[k],bin_maul_moment[k],color='#FFC125',height=w)
        plt.barh(bin_edges[k],bin_fault_moment[k],color='#EE2C2C',height=w)
        plt.barh(bin_edges[k],bin_moment[k],color='#3D59AB',height=w)
    elif bin_moment[k]>bin_maul_moment[k]:
        plt.barh(bin_edges[k],bin_fault_moment[k],color='#EE2C2C',height=w)
        plt.barh(bin_edges[k],bin_val_moment[k],color='#228B22',height=w)
        plt.barh(bin_edges[k],bin_maul_moment[k],color='#FFC125',height=w)
        plt.barh(bin_edges[k],bin_moment[k],color='#3D59AB',height=w)       
    else:
        plt.barh(bin_edges[k],bin_maul_moment[k],color='#FFC125',height=w)
        plt.barh(bin_edges[k],bin_val_moment[k],color='#228B22',height=w)
        plt.barh(bin_edges[k],bin_fault_moment[k],color='#EE2C2C',height=w)
        plt.barh(bin_edges[k],bin_moment[k],color='#3D59AB',height=w)
        
plt.ylim([bin_edges.min(),bin_edges.max()])
plt.xlim([1e17,1e22])
ax.set_xticks([1e17,1e19,1e21])
ax.yaxis.set_ticklabels([])
plt.xlabel('Cumulative Moment (N-m)')
plt.ylabel('Latitude')
plt.subplots_adjust(left=0.26,bottom=0.07,top=0.87,right=0.85)

plt.show()