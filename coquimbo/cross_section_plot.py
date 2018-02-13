'''
Make cross section plot with moment tensors, slip and slab
'''

from numpy import savetxt,arange,interp,linspace,array,c_,genfromtxt,expand_dims,zeros,cos,sin,deg2rad,where
import subprocess
from shlex import split
from matplotlib import pyplot as plt
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.imaging.mopad_wrapper import Beach
from mudpy import analysis,gmttools
from obspy.imaging.scripts import mopad
from matplotlib import cm as cmap
from matplotlib.ticker import AutoMinorLocator
import matplotlib
matplotlib.rcParams['font.size'] = 16


#P1
lonP1=[-74.1742,-70.208263]
latP1=[-30.2722,-30.660830]
#P2
lonP2=[-74.2646,-70.5483]
latP2=[-30.953,-31.4267]
#P3
lonP3=[-74.2458,-70.4235]
latP3=[-31.5494,-32.2538]
#ALL
Plon=[lonP1,lonP2,lonP3]
Plat=[latP1,latP2,latP3]

table='/Users/dmelgar/code/GMT/coquimbo/slice_table.txt'

#make colormap
#cm=gmttools.gmtColormap(u'/Users/dmelgar/code/python/cpt/color_linear.cpt')
cm=cmap.jet

xl=[-20,170]

#Start figure
fig, axarr = plt.subplots(len(Plon),1,figsize=(9,8)) 

for kP in range(len(Plon)):
    xi=linspace(Plon[kP][0],Plon[kP][1],1000)
    yi=interp(xi, array([Plon[kP][0],Plon[kP][1]]), array([Plat[kP][0],Plat[kP][1]]))
    savetxt(table,c_[xi,yi],fmt='%.6f\t%.6f')
    
    #run gmt slicing
    run='/Users/dmelgar/code/GMT/coquimbo/slice_slip.gmt'
    run=split(run)
    p=subprocess.Popen(run)
    p.communicate()
    
    #get the slab
    slab=genfromtxt('/Users/dmelgar/code/GMT/coquimbo/slab.txt')
    distance=zeros(len(slab))
    #Convert to distance
    for k in range(len(slab)):
        distance[k],az,baz=gps2DistAzimuth(slab[0,1],slab[0,0],slab[k,1],slab[k,0])
    distance=distance/1000
    
    
    #Get the topo/bathy
    topo=genfromtxt('/Users/dmelgar/code/GMT/coquimbo/topo.txt')
    distance_topo=zeros(len(topo))
    #Convert to distance
    for k in range(len(topo)):
        distance_topo[k],az,baz=gps2DistAzimuth(slab[0,1],slab[0,0],topo[k,1],topo[k,0])
        if az>baz:
            distance_topo[k]=-distance_topo[k]
    distance_topo=distance_topo/1000
    topo=topo[:,2]/1000
    
    #get the slip
    slip=genfromtxt('/Users/dmelgar/code/GMT/coquimbo/slip_profile.txt')
    
    #Get the moment tensors
    MTs=genfromtxt('/Users/dmelgar/code/GMT/coquimbo/MTs_slice.txt')
    #Convert to distance
    distance_MTs=zeros(len(MTs))
    for k in range(len(MTs)):
        distance_MTs[k],az,baz=gps2DistAzimuth(slab[0,1],slab[0,0],MTs[k,1],MTs[k,0])
        if az>baz:
            distance_MTs[k]=-distance_MTs[k]
    distance_MTs=distance_MTs/1000
    depth_MTs=-MTs[:,2]
    
    
    #Plot
    ax=axarr[kP]
    marker_color=expand_dims(slip[:,2],1)
    i=where(slab[:,2]>-60)[0]
    x=distance[i]
    y=slab[i,2]
    z=marker_color[i]
    im=ax.scatter(x,y,marker='o',lw=0,c=z,s=90,vmin=0,vmax=8,cmap=cm)
    # ADD MTs
    for k in range(len(MTs)):
        mt=MTs[k,3:9]
        #MT=analysis.MT(mt[0],mt[1],mt[2],mt[3],mt[4],mt[5],1,1,1)
        #MT.get_nodal_planes()
        MT=mopad.MomentTensor(mt,system='USE')
        np1,np2=MT.get_fps()
        #Thrust
        if (np1[2]>45 and np1[2]<135) or (np2[2]>45 and np2[2]<135):
            print k
            beach=Beach(np1, xy=(distance_MTs[k], depth_MTs[k]), width=6,linewidth=0.5,facecolor='r')
            print 'Thrust: '+str(np1)+' , '+str(np2)
        #Normal
        elif (np1[2]>-135 and np1[2]<-45) or (np2[2]>-135 and np2[2]<-45):
            print k
            beach=Beach(np1, xy=(distance_MTs[k], depth_MTs[k]), width=6,linewidth=0.5,facecolor='b')
            print 'Normal: '+str(np1)+' , '+str(np2)
        #Oblique
        else:
            print k
            beach=Beach(np1, xy=(distance_MTs[k], depth_MTs[k]), width=6,linewidth=0.5,facecolor='g')
            print 'Oblique: '+str(np1)+' , '+str(np2)
        ax.add_collection(beach) 
    #plt.colorbar()
    ax.plot(distance_topo,topo,lw=2)
    ax.yaxis.tick_left()
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_label_position("right")
    ax.set_ylim([-60,10])
    ax.set_yticklabels(['','50','40','30','20','10','0',''])
    ax.set_ylabel('Depth (km)',fontsize=16)
    ax.set_xlim(xl)
    if kP!=2:
        ax.xaxis.set_ticklabels([])
    if kP==0:
        ax.annotate("a",xy=(-15,2),annotation_clip=False)
        ax.annotate("a'",xy=(160,2),annotation_clip=False)
    if kP==1:
        ax.annotate("b",xy=(-15,2),annotation_clip=False)
        ax.annotate("b'",xy=(160,2),annotation_clip=False)
    if kP==2:
        ax.annotate("c",xy=(-15,2),annotation_clip=False)
        ax.annotate("c'",xy=(160,2),annotation_clip=False)
        ax.set_xlabel('Distance from the trench (km)',fontsize=16)
    
    minorLocatorx = AutoMinorLocator()
    ax.xaxis.set_minor_locator(minorLocatorx)
    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=3) 
    minorLocatory = AutoMinorLocator()   
    ax.yaxis.set_minor_locator(minorLocatory)
    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=3)  
    
    #ax.title('Lower Hemisphere Projection')
    
    
    #plt.figure(figsize=(14,4))
    #ax=plt.gca()
    #marker_color=expand_dims(slip[:,2],1)
    #plt.scatter(distance,slab[:,2],marker='o',lw=0,c=marker_color,s=60,vmin=0,vmax=10,cmap=cm)
    ##Rotation matrix
    #th=deg2rad(-90)
    ##Ry=array([[cos(th),0,sin(th)],[0,1,0],[-sin(th),0,cos(th)]])
    ##Rx=array([[1,0,0],[0,cos(th),-sin(th)],[0,sin(th),cos(th)]])
    #Rz=array([[cos(th),-sin(th),0],[sin(th),cos(th),0],[0,0,1]])
    ## ADD MTs
    #for k in range(len(MTs)):
    #    mt=MTs[k,3:9]
    #    mt=array([[mt[0],mt[3],mt[4]],[mt[3],mt[1],mt[5]],[mt[4],mt[5],mt[2]]])
    #    MT=mopad.MomentTensor([mt[0,0],mt[1,1],mt[2,2],mt[0,1],mt[0,2],mt[1,2]],system='USE')
    #    np1,np2=MT.get_fps()
    #    #Now make the roated one for plotting
    #    mt=Rz.dot(mt)
    #    MTplot=mopad.MomentTensor([mt[0,0],mt[1,1],mt[2,2],mt[0,1],mt[0,2],mt[1,2]],system='USE')
    #    np1_plot,np2_plot=MTplot.get_fps()
    #    #Thrust
    #    if (np1[2]>45 and np1[2]<135) or (np2[2]>45 and np2[2]<135):
    #        beach=Beach(np1_plot, xy=(distance_MTs[k], depth_MTs[k]), width=5,linewidth=0.5,facecolor='r')
    #    #Normal
    #    elif (np1[2]>-135 and np1[2]<-45) or (np2[2]>-135 and np2[2]<-45):
    #        beach=Beach(np1_plot, xy=(distance_MTs[k], depth_MTs[k]), width=5,linewidth=0.5,facecolor='b')
    #    #Oblique
    #    else:
    #        beach=Beach(np1_plot, xy=(distance_MTs[k], depth_MTs[k]), width=5,linewidth=0.5,facecolor='g')
    #    ax.add_collection(beach) 
    #plt.colorbar()
    #plt.plot(distance_topo,topo,lw=2)
    #plt.xlabel('Distance from the trench (km)',fontsize=16)
    #plt.ylabel('Depth (km)',fontsize=16)
    #plt.ylim([-60,5])
    #plt.xlim([-20,200])
    #plt.subplots_adjust(top=0.8,bottom=0.2)
    #plt.title('Right-hand Side Hemisphere Projection')

plt.subplots_adjust(top=0.9,bottom=0.1,right=0.85,left=0.15,hspace=0.02)
#Make colrobar
cbar_ax = fig.add_axes([0.15, 0.91, 0.7, 0.02])
cb=fig.colorbar(im, cax=cbar_ax,orientation='horizontal')
cb.set_label('Slip (m)', labelpad=-50)
for t in cbar_ax.xaxis.get_major_ticks(): 
    t.tick1On = True 
    t.tick2On = True 
    t.label1On = False 
    t.label2On = True 

plt.show()