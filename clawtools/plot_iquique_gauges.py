from matplotlib import pyplot as plt
from numpy import genfromtxt,unique,zeros,where,array,nan,c_
from obspy import read

gauges_file=u'/Users/dmelgar/Tsunamis/iquique_42_kinematic/_output/fort.gauge'
gauges_dict=u'/Users/dmelgar/Slip_inv/iquique_42/data/station_info/gauges.dict'
pathdata=u'/Users/dmelgar/Iquique2014/tsunami/gauge_data/'


clawname=genfromtxt(gauges_dict,usecols=0)
gaugename=genfromtxt(gauges_dict,usecols=3,dtype='S')

all_syn=genfromtxt(gauges_file)


fig, axarr = plt.subplots(4, 1)  
for k in range(len(clawname)):
    ax=axarr[k]
    gauge=read(pathdata+gaugename[k]+'.sac')
    isyn=where(all_syn[:,0]==clawname[k])[0]
    tsyn=all_syn[isyn,2]
    syn=all_syn[isyn,-1]
    ax.plot(gauge[0].times()/60,gauge[0].data,'k')
    ax.plot(tsyn/60,syn,'r')
    ax.set_ylabel(gaugename[k]+' (m)')
    if k==3:
        plt.xlabel('Minutes after origin time')
    ax.grid()
plt.subplots_adjust(left=0.2, bottom=0.05, right=0.8, top=0.95)
    
plt.show()
