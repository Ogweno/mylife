from numpy import genfromtxt,arange,diff,where,log10,isnan
from scipy.stats import binned_statistic
from matplotlib import pyplot as plt

pop_density=genfromtxt('/Users/dmelgar/code/GMT/Mexico_GMPEs/log_people.xyz')
slope=genfromtxt('/Users/dmelgar/DEMs/srtm15/central_mexico_slope_at_people.xyz')

#clean up
i=where(pop_density[:,2]<0)[0]
pop_density[:,2][i]=0
i=where(pop_density[:,1]<14.5293)[0]
pop_density[:,2][i]=0
i=where(isnan(slope[:,2])==True)[0]
slope[i,2]=1e3
i=where(slope[:,2]<=0)[0]
slope[i,2]=1e-6
i=where(slope[:,2]<1e-6)[0]
slope[i,2]=1e-6





bins=arange(1,4.6,0.05)
stat,edges,numbers=binned_statistic(pop_density[:,2],slope[:,2],bins=bins,statistic='median')
count,edges,numbers=binned_statistic(pop_density[:,2],slope[:,2],bins=bins,statistic='count')
centers=diff(edges)+edges[0:-1]
total=(0.45**2)*count*10**(centers)



plt.figure(figsize=(5,4))
plt.bar(edges[0:-1],stat,width=0.05,facecolor='#B0C4DE')

ax=plt.gca()
labels=['',r'$10^1$','',r'$10^2$','',r'$10^3$','',r'$10^4$','']
ax.xaxis.set_ticklabels(labels,fontsize=14)

plt.xlabel(r'Population density (people/km$^2$)',fontsize=14)
plt.ylabel('Median slope (%)',fontsize=14)
plt.yticks(fontsize=14)
plt.grid()
ax.set_axisbelow(True)
plt.xlim([0.94,4.6])
plt.subplots_adjust(bottom=0.2,left=0.2)







##Calculate people per slope bin
bin_width=0.2
total_people_in_pixel=(10**pop_density[:,2])*(0.45**2)
slope_bins=arange(0,20,bin_width)
stat,edges,numbers=binned_statistic(slope[:,2],total_people_in_pixel,bins=slope_bins,statistic='sum')

plt.figure(figsize=(5,4))
plt.bar(edges[0:-1],stat/1e6,width=bin_width,facecolor='#B0C4DE')

ax=plt.gca()
labels=['','0','2','4','6','8','10']
ax.xaxis.set_ticklabels(labels,fontsize=14)

plt.xlabel(r'Slope (%)',fontsize=14)
plt.ylabel('People (millions)',fontsize=14)
plt.yticks(fontsize=14)
plt.grid()
ax.set_axisbelow(True)
plt.xlim([-0.2,10])
plt.subplots_adjust(bottom=0.2,left=0.2)




plt.show()