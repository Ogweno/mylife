from numpy import array,zeros,genfromtxt,arange,sqrt
from matplotlib import pyplot as plt
import matplotlib

#Make the colormaps for plotting later on
cdict = {'red': ((0., 1, 1),
                 (0.10, 1, 1),
                 (0.20, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.10, 1, 1),
                   (0.20, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.15, 1, 1),
                  (0.20, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}
whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)



#Read in displacememnts for the different dip angles
d65=genfromtxt('/Users/dmelgar/code/MFILES/Coulomb3/output_files/lefkada_synth65.displ')
d70=genfromtxt('/Users/dmelgar/code/MFILES/Coulomb3/output_files/lefkada_synth70.displ')
d75=genfromtxt('/Users/dmelgar/code/MFILES/Coulomb3/output_files/lefkada_synth75.displ')

#Define the stallite look direction UNIT vector (e,n,up)
look=array([-0.562,-0.105,0.820])

#Project displacements along LOS vector
los65=zeros(len(d65))
los70=zeros(len(d70))
los75=zeros(len(d75))
for k in range(len(d65)):
    los65[k]=(d65[k,3]*look[0]+d65[k,4]*look[1]+d65[k,5]*look[2])
    los70[k]=(d70[k,3]*look[0]+d70[k,4]*look[1]+d70[k,5]*look[2])
    los75[k]=(d75[k,3]*look[0]+d75[k,4]*look[1]+d75[k,5]*look[2])
 
    
    
    
    
    
#Make the plots    
i=arange(0,len(los65),78) 
fig=plt.figure(figsize=(7,7))

#plot quivers and colors first
ax=plt.subplot(231)
ax.scatter(d65[:,0],d65[:,1],c=sqrt(d65[:,3]**2+d65[:,4]**2),cmap=whitejet,vmin=0,vmax=0.4,lw=0)
ax.quiver(d65[i,0],d65[i,1],d65[i,3],d65[i,4],width=0.005,headwidth=4,scale=3)
ax.plot([-5.13,5.13],[-14.095,14.095],c='#101010',lw=3)
ax.set_xlim([-15,15])
ax.set_ylim([-25,25])
ax.xaxis.set_ticklabels([])
ax.set_ylabel('Y (km)')
ax.set_title(r'$\delta = 65^\circ$')

ax=plt.subplot(232)
ax.scatter(d70[:,0],d70[:,1],c=sqrt(d70[:,3]**2+d70[:,4]**2),cmap=whitejet,vmin=0,vmax=0.4,lw=0)
ax.quiver(d70[i,0],d70[i,1],d70[i,3],d70[i,4],width=0.005,headwidth=4,scale=3)
ax.plot([-5.13,5.13],[-14.095,14.095],c='#101010',lw=3)
ax.set_xlim([-15,15])
ax.set_ylim([-25,25])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.set_title(r'$\delta = 70^\circ$')

ax=plt.subplot(233)
cf0=ax.scatter(d75[:,0],d75[:,1],c=sqrt(d75[:,3]**2+d75[:,4]**2),cmap=whitejet,vmin=0,vmax=0.4,lw=0)
ax.quiver(d75[i,0],d75[i,1],d75[i,3],d75[i,4],width=0.005,headwidth=4,scale=3)
ax.plot([-5.13,5.13],[-14.095,14.095],c='#101010',lw=3)
ax.set_xlim([-15,15])
ax.set_ylim([-25,25])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.set_title(r'$\delta = 75^\circ$')

ax=plt.subplot(234)
ax.scatter(d65[:,0],d65[:,1],c=los65,cmap=plt.cm.seismic,vmin=-0.15,vmax=0.15,lw=0)
ax.plot([-5.13,5.13],[-14.095,14.095],c='#101010',lw=3)
ax.set_xlim([-15,15])
ax.set_ylim([-25,25])
ax.set_ylabel('Y (km)')
ax.xaxis.set_ticklabels(['','-10','','0','','10',''])

ax=plt.subplot(235)
ax.scatter(d70[:,0],d70[:,1],c=los70,cmap=plt.cm.seismic,vmin=-0.15,vmax=0.15,lw=0)
ax.plot([-5.13,5.13],[-14.095,14.095],c='#101010',lw=3)
ax.set_xlim([-15,15])
ax.set_ylim([-25,25])
ax.yaxis.set_ticklabels([])
ax.set_xlabel('X (km)')
ax.xaxis.set_ticklabels(['','-10','','0','','10',''])


ax=plt.subplot(236)
cf=ax.scatter(d75[:,0],d75[:,1],c=los75,cmap=plt.cm.seismic,vmin=-0.15,vmax=0.15,lw=0)
ax.plot([-5.13,5.13],[-14.095,14.095],c='#101010',lw=3)
ax.set_xlim([-15,15])
ax.set_ylim([-25,25])
ax.yaxis.set_ticklabels([])
ax.xaxis.set_ticklabels(['','-10','','0','','10',''])

#Colorbars
position=fig.add_axes([0.85,0.53,0.02,0.35])  ## the parameters are the specified position you set 
fig.colorbar(cf0,cax=position,label='Horizontal (m)',ticks=[0,0.1,0.2,0.3,0.4])
position=fig.add_axes([0.85,0.12,0.02,0.35])  ## the parameters are the specified position you set 
fig.colorbar(cf,cax=position,label='LOS (m)',ticks=[-0.15,-0.10,-0.05,0,0.05,0.1,0.15])


plt.subplots_adjust(wspace=0.05,hspace=0.05,right=0.84)
plt.show()