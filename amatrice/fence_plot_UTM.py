from numpy import genfromtxt,array,zeros,where,r_,arange,mean
from matplotlib import pyplot as plt
from pyproj import Proj
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator
from mudpy import gmttools
from mudpy.viewFQ import get_subfault_corners
from pyproj import Proj

faults_to_plot=5# 1 to 5main,synth,pian,anti,detach
pian=False
#Run time parameters
home='/Users/dmelgar/Slip_inv/'

project_name='Amatrice_M6.6'
#project_name='Amatrice_M6.1_August'
#project_name='Amatrice_M6.1_October'

thresh=0

#run_name='gps_insar_vettore_main_synth_pian_anti_detach'
#run_number='0022'
#run_name='gps_insar_vettore_main'
#run_number='0017'
#run_name='gps_insar_vettore_main_pian'
#run_number='0012'
run_name='jacknife.std'
run_number='0000'

meshfile1='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_main.mshout'
meshfile2='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_synth.mshout'
meshfile3='/Users/dmelgar/Amatrice2016/3D_fault/_complex/pian_grande.mshout'
meshfile4='/Users/dmelgar/Amatrice2016/3D_fault/_complex/antithetic_listric.mshout'

#meshfile1='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_main.mshout'
#meshfile2='/Users/dmelgar/Amatrice2016/3D_fault/_complex/pian_grande.mshout'
#meshfile3='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_synth.mshout'
#meshfile4='/Users/dmelgar/Amatrice2016/3D_fault/_complex/antithetic_listric.mshout'

plt.rcParams['font.size'] = 10
af=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records.txt',usecols=[0,1,2,3,4,5])
af=array([[-100,-100,-100,-100,-100]])
maxslip=None
UTM_zone='33T'
fudge=0.01
fake_hypo=[0.1,0.1],
borderwidth=0.1
figsize=(12,3)
xtick=20
ytick=20
ztick=2
inverse_model=False
hypocenter=None
strike=155
xl=[0,30]
xorigin=325
yorigin=4710
yl=[0,50]
zl=[0,10]

#Proj stuff
p = Proj(proj='utm',zone=UTM_zone,ellps='WGS84')


fault_name=home+project_name+'/output/inverse_models/models/%s.%s.inv' % (run_name,run_number)
gmttools.make_total_model(fault_name,thresh=0)
fault=genfromtxt(home+project_name+'/output/inverse_models/models/%s.%s.inv.total' % (run_name,run_number))
reference=genfromtxt('/Users/dmelgar/Slip_inv/Amatrice_M6.6/output/inverse_models/models/gps_insar_vettore_main_synth_pian_anti_detach.0018.inv')
#Parse log file for hypocenter
log_file=home+project_name+'/output/inverse_models/models/%s.%s.log' % (run_name,run_number)
f=open(log_file,'r')
loop_go=True
while loop_go:
    line=f.readline()  
    if 'Mw' in line:
        Mw=float(line.split(':')[-1].split(' ')[-1])   
        loop_go=False
f.close() 
    

def get_corners(meshfile):
    #get subfault corners
    corners=genfromtxt(meshfile,usecols=range(4,13))
    i=where(corners[:,0]>360)[0]
    corners[i,0]=corners[i,0]-360
    i=where(corners[:,3]>360)[0]
    corners[i,3]=corners[i,3]-360
    i=where(corners[:,6]>360)[0]
    corners[i,6]=corners[i,6]-360
    corners[:,2]=-corners[:,2]
    corners[:,5]=-corners[:,5]
    corners[:,8]=-corners[:,8]
    #Convert to UTM kilometers
    x,y=p(corners[:,0],corners[:,1])
    corners[:,0]=x/1000-xorigin
    corners[:,1]=y/1000-yorigin
    x,y=p(corners[:,3],corners[:,4])
    corners[:,3]=x/1000-xorigin
    corners[:,4]=y/1000-yorigin
    x,y=p(corners[:,6],corners[:,7])
    corners[:,6]=x/1000-xorigin
    corners[:,7]=y/1000-yorigin
    return corners
                
corners1=get_corners(meshfile1)  
corners2=get_corners(meshfile2) 
corners3=get_corners(meshfile3) 
corners4=get_corners(meshfile4)  
if faults_to_plot==1:
    corners=corners1
elif faults_to_plot==2:    
    corners=r_[corners1,corners2]
elif faults_to_plot==3:    
    corners=r_[corners1,corners2,corners3]
elif faults_to_plot>=4:    
    corners=r_[corners1,corners2,corners3,corners4]

Nfaults=len(corners)
    
    
#Normalized slip
slip=(fault[:,8]**2+fault[:,9]**2)**0.5
i=where(slip<thresh)[0]
slip[i]=0
total_slip=slip

#Saturate to maxslip
if maxslip!=None:
    imax=where(slip>maxslip)[0]
    slip[imax]=maxslip
#normalize
norm_slip=slip/slip.max()

#get detachment corners
if faults_to_plot==5:
    fault_detach=fault[len(corners1)+len(corners2)+len(corners3)+len(corners4):-1]
else:
    fault_detach=reference[len(corners1)+len(corners2)+len(corners3)+len(corners4):-1]
nrows=4
fault_detach=fault_detach[0:-nrows*17+1]
fault_detach=fault_detach[2*17:,:]
slip_detach=norm_slip[len(corners1)+len(corners2)+len(corners3)+len(corners4):-1]
corners_detach=get_subfault_corners(fault_detach)
x,y=p(corners_detach[:,0],corners_detach[:,1])
corners_detach[:,0]=x/1000-xorigin
corners_detach[:,1]=y/1000-yorigin
x,y=p(corners_detach[:,3],corners_detach[:,4])
corners_detach[:,3]=x/1000-xorigin
corners_detach[:,4]=y/1000-yorigin
x,y=p(corners_detach[:,6],corners_detach[:,7])
corners_detach[:,6]=x/1000-xorigin
corners_detach[:,7]=y/1000-yorigin
x,y=p(corners_detach[:,9],corners_detach[:,10])
corners_detach[:,9]=x/1000-xorigin
corners_detach[:,10]=y/1000-yorigin

#Get colormaps
#slip_colormap = colors.LinearSegmentedColormap('slip',whitejet_dict,256)
slip_colormap = plt.cm.gist_heat_r


def axisEqual3D(ax):
    extents = array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


#Plot init, axes positions etc
def make_subplot(ax1,sub_corners,sub_slip,af,azimuth,elevation,title,colorbar,triangles=True,axison=False,xl=[],yl=[],axisequal=False):
    
    ax1.set_xlim(xl)
    ax1.set_ylim(yl)
    ax1.set_zlim(zl)
    #Fenagle the axis ticks
    xmajorLocator = MultipleLocator(xtick)
    ymajorLocator = MultipleLocator(ytick)
    zmajorLocator = MultipleLocator(ztick)
    ax1.xaxis.set_major_locator(xmajorLocator)
    ax1.yaxis.set_major_locator(ymajorLocator)
    ax1.zaxis.set_major_locator(zmajorLocator)
    ax1.invert_zaxis()
    ax1.view_init(elev=elevation, azim=azimuth)
    
    #Make one patch per subfaultr
    for ksub in range(len(sub_corners)):
        if triangles:
            vertices=[[tuple(sub_corners[ksub,0:3]),tuple(sub_corners[ksub,3:6]),tuple(sub_corners[ksub,6:9])]]
        else:
            vertices=[[tuple(sub_corners[ksub,0:3]),tuple(sub_corners[ksub,3:6]),tuple(sub_corners[ksub,6:9]),tuple(sub_corners[ksub,9:12])]]
        subfault=Poly3DCollection(vertices, linewidths=borderwidth)
        #subfault.set_color(pqlx(norm_slip[ksub]))
        subfault.set_color(slip_colormap(sub_slip[ksub]))
        subfault.set_linewidth(borderwidth)
        subfault.set_edgecolor('#505050')
        ax1.add_collection3d(subfault)
    
    #Add aftershocks
    #ax1.scatter(af[:,3],af[:,2],af[:,4],s=5,lw=0)
    
    if colorbar==True:
        #Dummy mapable for colorbar
        s=plt.scatter(zeros(len(total_slip)),zeros(len(total_slip)),c=total_slip,cmap=slip_colormap,s=0.00001,lw=0)
        
        #Mke colorbar
        cb=plt.colorbar(s,shrink=0.9,pad=0.1)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator=tick_locator
        cb.update_ticks()
        cb.set_label('Slip (m)')
    
    #Labels n' stuff
    ax1.set_xlabel('\n\nEast (km)')
    ax1.set_ylabel('\n\nNorth (km)')
    ax1.set_zlabel('z(km)',rotation=90)
    ax1.set_title(title)
    if axison==False:
        ax1.set_axis_off()
    
def move_corners(corners,corner_type=3):
    
    if corner_type==3:
    
        delta_x=r_[corners[:,0],corners[:,3],corners[:,6]].min()
        delta_y=r_[corners[:,1],corners[:,4],corners[:,7]].min()
        corners[:,0]=corners[:,0]-delta_x
        corners[:,3]=corners[:,3]-delta_x
        corners[:,6]=corners[:,6]-delta_x
        corners[:,1]=corners[:,1]-delta_y
        corners[:,4]=corners[:,4]-delta_y
        corners[:,7]=corners[:,7]-delta_y
    
    else:
    
        delta_x=r_[corners[:,0],corners[:,3],corners[:,6],corners[:,9]].min()
        delta_y=r_[corners[:,1],corners[:,4],corners[:,7],corners[:,10]].min()
        corners[:,0]=corners[:,0]-delta_x
        corners[:,3]=corners[:,3]-delta_x
        corners[:,6]=corners[:,6]-delta_x
        corners[:,9]=corners[:,9]-delta_x
        corners[:,1]=corners[:,1]-delta_y
        corners[:,4]=corners[:,4]-delta_y
        corners[:,7]=corners[:,7]-delta_y      
        corners[:,10]=corners[:,10]-delta_y
    
    return corners

 
            

#Determine whch indices for which fault
i_vettore_main=arange(0,len(corners1))
i_vettore_synth=arange(len(corners1),len(corners1)+len(corners2))
if pian==False:
    i_pian=arange(len(corners1)+len(corners2),len(corners1)+len(corners2)+len(corners3))
else:
    i_pian=arange(len(corners1),len(corners1)+len(corners3))
i_anti=arange(len(corners1)+len(corners2)+len(corners3),len(corners1)+len(corners2)+len(corners3)+len(corners4))

fig=plt.figure(figsize=figsize)

#Focus on Vettore main
xl=[5,30]
yl=[7,53]
ax = fig.add_subplot(162, projection='3d')
plot_corners=move_corners(corners[i_vettore_main])
make_subplot(ax,plot_corners,norm_slip[i_vettore_main],af,azimuth=-100,elevation=57,title='MVFM',colorbar=False,xl=xl,yl=yl)

#Focus on vettore synth
xl=[6,24]
yl=[11,50]
ax = fig.add_subplot(163, projection='3d')
if faults_to_plot>=2:
    if pian==False:
        plot_corners=move_corners(corners[i_vettore_synth])
        make_subplot(ax,plot_corners,norm_slip[i_vettore_synth],af,azimuth=-100,elevation=57,title='MVFS',colorbar=False,xl=xl,yl=yl)
else:
    plot_corners=move_corners(corners2)
    make_subplot(ax,plot_corners,zeros(len(plot_corners)),af,azimuth=-100,elevation=57,title='MVFS',colorbar=False,xl=xl,yl=yl)        

#Focus on pian grande
xl=[0,18]
yl=[12,45]
ax = fig.add_subplot(164, projection='3d')
if faults_to_plot>=3:
    plot_corners=move_corners(corners3)
    make_subplot(ax,plot_corners,norm_slip[i_pian],af,azimuth=-10,elevation=46,title='CAF',colorbar=False,xl=xl,yl=yl)
else:
    plot_corners=move_corners(corners3)
    make_subplot(ax,plot_corners,zeros(len(plot_corners)),af,azimuth=-10,elevation=77,title='CAF',colorbar=False,triangles=True,xl=xl,yl=yl,axison=False)
    
#Focus on antithetic
xl=[0,30]
yl=[10,50]
ax = fig.add_subplot(165, projection='3d')
if faults_to_plot>=4:
    plot_corners=move_corners(corners[i_anti])
    make_subplot(ax,plot_corners,norm_slip[i_anti],af,azimuth=-10,elevation=46,title='SAF',colorbar=False,xl=xl,yl=yl)
else:
    plot_corners=move_corners(corners4)
    make_subplot(ax,plot_corners,zeros(len(plot_corners)),af,azimuth=-10,elevation=77,title='SAF',colorbar=False,triangles=True,xl=xl,yl=yl,axison=False)

#Focus on detachment
xl=[7,45]
yl=[13,58]
ax = fig.add_subplot(166, projection='3d')
plot_corners=move_corners(corners_detach,corner_type=4)
if faults_to_plot==5:
    make_subplot(ax,plot_corners,slip_detach,af,azimuth=-10,elevation=77,title='AF',colorbar=True,triangles=False,xl=xl,yl=yl,axison=False)
else:
    make_subplot(ax,plot_corners,zeros(len(plot_corners)),af,azimuth=-10,elevation=77,title='AF',colorbar=True,triangles=False,xl=xl,yl=yl,axison=False)



#Plot all faults
#xl=[0,50]
#yl=[-30,65]
#ax = fig.add_subplot(161, projection='3d')
#if faults_to_plot==5:
#    make_subplot(ax,corners_detach,slip_detach,af,azimuth=-60,elevation=51,title='',colorbar=False,triangles=False,axison=True,xl=xl,yl=yl)
#if pian==False:
#    make_subplot(ax,corners,norm_slip,af,azimuth=-60,elevation=51,title='All',colorbar=False,axison=True,axisequal=True,xl=xl,yl=yl)
#    print 'a'


plt.subplots_adjust(hspace=0.2,wspace=0.2,left=0,right=0.95,top=0.9,bottom=0.22)

plt.show()