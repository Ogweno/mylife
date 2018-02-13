from numpy import genfromtxt,zeros
from matplotlib import pyplot as plt

#file_list=['/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_20km.txt',
#    '/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_30km.txt',
#    '/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_40km.txt',
#    '/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_50km.txt',
#    '/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_75km.txt',
#    '/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_100km.txt']
    
file_list=['/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_20km.txt',
    '/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_30km.txt',
    '/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_40km.txt',
    '/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_50km.txt',
    '/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_75km.txt',
    '/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_100km.txt']
    
spacing=[20,30,40,50,75,100]
    
sta=zeros(len(file_list))
for k in range(len(file_list)):
    s=genfromtxt(file_list[k],usecols=2)
    sta[k]=s.mean()
    
plt.figure()
plt.plot(spacing,sta,lw=1.5,c='k')
plt.plot([0,110],[6,6],'--',lw=1.0)
plt.scatter(spacing,sta,s=50,marker='s',color='r')
plt.xlabel('Station spacing (km)')
plt.ylabel('# of stations sensed event')
plt.ylim([0,30])
plt.xlim([18,105])
plt.show()
