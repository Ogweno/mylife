from numpy import genfromtxt,log10,histogram,diff,linspace
from matplotlib import pyplot as plt
from scipy.stats import binned_statistic

stagps='RCSD'
catalog=genfromtxt('/Users/dmelgar/Valparaiso2017/catalog/catalog.txt')
lonlat=genfromtxt('/Users/dmelgar/Valparaiso2017/GPS/gps.sta',usecols=[1,2])
sta=genfromtxt('/Users/dmelgar/Valparaiso2017/GPS/gps.sta',usecols=0,dtype='S')
t0=4*30+22+22/24.+46/3600.
days=catalog[:,1]*30+catalog[:,2]+catalog[:,3]/24.+catalog[:,4]/3600.
days=days-t0

gps=genfromtxt('/Users/dmelgar/Valparaiso2017/GPS/'+stagps+'.txt')
tgps=gps[:,2]-112-22/24.-46/3600. #Aprill 22nd is dat 112
disp=gps[:,3]

#Map
plt.figure(figsize=(6,18))
plt.scatter(catalog[:,7],catalog[:,6],marker='o',s=3**catalog[:,11],lw=1)
plt.scatter(lonlat[:,0],lonlat[:,1],marker='^',s=90,c='r')
#plt.axis('equal')
for k in range(len(sta)):
    if sta[k]=='TRPD':
        plt.annotate(s=sta[k],xy=(lonlat[k,0]+0.07,lonlat[k,1]+0.05),fontsize=14)
    elif sta[k]=='VALN':
        plt.annotate(s=sta[k],xy=(lonlat[k,0]+0.07,lonlat[k,1]-0.01),fontsize=14)
    else:
        plt.annotate(s=sta[k],xy=(lonlat[k,0]+0.07,lonlat[k,1]-0.07),fontsize=14)


#Statistics
plt.figure(figsize=(20,5))
plt.subplot(121)
markerline, stemlines, baseline=plt.stem(days,catalog[:,11])
plt.setp(stemlines, 'color', 'k',lw=0.1)
plt.setp(markerline,'color','r','markersize',5)
plt.ylim([2,7.2])
plt.xlim([-0.5,14])
plt.xlabel('Days since start of swarm')
plt.ylabel('Magnitude')

ax=plt.subplot(122)
nbins=15
ax.set_yscale('log')
n, bin_edges = histogram(catalog[:,11], nbins)
bin_center = diff(bin_edges)/2.+bin_edges[0:-1]
ax.scatter(bin_center,n,c='r',s=50)
ax.hist(catalog[:,11],nbins,histtype='step')

mag=linspace(0,7.2)
freq=10**(-mag+5.5)
plt.plot(mag,freq,'--',lw=2,c='k')
ax.set_xlabel('Magnitude')
ax.set_ylabel('No. of events')
ax.set_ylim([0.5,500])
ax.set_xlim([2,8])


#plot gps vs moment
moment=10**(1.5*catalog[:,11]+9.1)
stat,bin_edges,bin_num=binned_statistic(days, moment, statistic='sum', bins=20)
fig=plt.figure(figsize=(20,3))
ax1=fig.add_subplot(111)
ax1.bar(bin_edges[0:-1],log10(stat),width=0.6,color='#FF7F50')
ax1.set_ylim([15,21])
ax1.set_ylabel('Moment (Nm)')
ax1.set_xlabel('Days since start of swarm')
ax1.set_xlim([-100,20])

ax2 = ax1.twinx()
ax2.scatter(tgps,disp)
ax2.plot(tgps,disp)
ax2.set_ylabel('East (mm)')
ax2.plot([0,0],[-90,90],'--',c='k',lw=2)
ax2.set_ylim([-60,10])
plt.title('Station '+stagps)
ax2.set_xlim([-100,20])

plt.subplots_adjust(bottom=0.15)




plt.show()