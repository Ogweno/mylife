'''
Compute travel times for all PGD earthquakes
'''

from glob import glob
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import locations2degrees
from numpy import genfromtxt,float64,zeros,where,arange
from matplotlib import pyplot as plt
from string import rjust

make_travel=False
plot_pgd_dist=True
plot_pgd_dist_all=True

path='/Users/dmelgar/PGD/'
earthquakes=glob(path+'/station_info/*')
if make_travel:
    #initalize
    for k in range(len(earthquakes)):
        event=earthquakes[k].split('/')[-1].split('.')[0]
        print event
        #Get hypocenter
        hypo=genfromtxt(path+'event_info/'+event+'.hypo')
        #Read station data
        stanames=genfromtxt(earthquakes[k],usecols=0,dtype='S')
        station_coords=genfromtxt(earthquakes[k],usecols=[1,2])
        #initalize
        tp=zeros(len(stanames))
        delta=zeros(len(stanames))
        for j in range(len(stanames)):
            #Get event-station distance in degrees
            delta[j]=locations2degrees(station_coords[j,1],station_coords[j,0],hypo[2],hypo[1])
            #Get p-time to site
            tt=getTravelTimes(delta[j],hypo[3])
            tp[j]=float64(tt[0]['time'])
        #Write to file
        f=open(path+'travel_times/'+event+'.tt','w')
        for j in range(len(tp)):
            line='%s\t%8.2f\n' %(stanames[j],tp[j])
            f.write(line)
        f.close()
    
    
#Plot PGD computation as a function fo distance
if plot_pgd_dist:
        for k in [0,1,2,3,4,5,6]:#range(len(earthquakes)):
            event=earthquakes[k].split('/')[-1].split('.')[0]
            print event
            #Get stations and p-times
            stations=genfromtxt(path+'PGD/'+event+'.txt',usecols=0,dtype='S')
            tpgd=genfromtxt(path+'PGD/'+event+'.txt',usecols=5)
            distance=genfromtxt(path+'PGD/'+event+'.txt',usecols=1)
            sta_compare=genfromtxt(path+'travel_times/'+event+'.tt',usecols=0,dtype='S')
            tp=genfromtxt(path+'travel_times/'+event+'.tt',usecols=1)
            tfinal=zeros(len(tpgd))
            if event=='tohoku':
                for j in range(len(sta_compare)):
                    sta_compare[j]=rjust(sta_compare[j],4,'0')
            if event=='napa':
                for j in range(len(stations)):
                    i=where(sta_compare==stations[j].upper())[0]
                    tfinal[j]=tpgd[j]#-tp[i]
            else:
                for j in range(len(stations)):
                    i=where(sta_compare==stations[j])[0]
                    tfinal[j]=tpgd[j]#-tp[i]
            plt.figure()
            plt.plot(distance,tfinal,lw=0,marker='o')
            plt.grid()
            plt.xlabel('Hypocentral distance (km)')
            plt.ylabel('Time to PGD (s)')
            
            #Plot reference lines
            for j in range(5):
                d=arange(0,distance.max()+10,1)
                t=d/(1*j)
                plt.plot(d,t,'--',c='k')
            plt.xlim([0,distance.max()])
            plt.ylim([0,tfinal.max()])
            plt.title(event)
            plt.savefig(path+'plots/'+event+'.time2pgd.horiz.noP.png')
        plt.close("all")
        
if plot_pgd_dist_all:
        plt.figure()
        for k in [0,1,2,3,4,5,6]:#range(len(earthquakes)):
            event=earthquakes[k].split('/')[-1].split('.')[0]
            print event
            #Get stations and p-times
            stations=genfromtxt(path+'PGD/'+event+'.txt',usecols=0,dtype='S')
            tpgd=genfromtxt(path+'PGD/'+event+'.txt',usecols=5)
            distance=genfromtxt(path+'PGD/'+event+'.txt',usecols=1)
            sta_compare=genfromtxt(path+'travel_times/'+event+'.tt',usecols=0,dtype='S')
            tp=genfromtxt(path+'travel_times/'+event+'.tt',usecols=1)
            tfinal=zeros(len(tpgd))
            if event=='tohoku':
                for j in range(len(sta_compare)):
                    sta_compare[j]=rjust(sta_compare[j],4,'0')
            if event=='napa':
                for j in range(len(stations)):
                    i=where(sta_compare==stations[j].upper())[0]
                    tfinal[j]=tpgd[j]#-tp[i]
            else:
                for j in range(len(stations)):
                    i=where(sta_compare==stations[j])[0]
                    tfinal[j]=tpgd[j]#-tp[i]
            plt.plot(distance,tfinal,lw=0,marker='s',markersize=5)
            plt.grid()
            plt.xlabel('Hypocentral distance (km)')
            plt.ylabel('Time to PGD (s)')
            plt.title('All events')
        plt.legend(['M7.2 El Mayor','M8.2 Iquique','M8.8 Maule','M7.8 Mentawai','M6.1 Napa','M6.0 Parkfield','M9.0 Tohoku'],loc=4,prop={'size':12})
        #Plot reference lines
        for j in range(5):
            d=arange(0,1010,1)
            t=d/(1*j)
            plt.plot(d,t,'--',c='k')
        plt.xlim([0,1000])
        plt.ylim([0,350])
        plt.savefig(path+'plots/allevents.time2pgd.horiz.noP.png')
        plt.close("all")

    

