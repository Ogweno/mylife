from numpy import genfromtxt,array,zeros,where,r_,arange
from matplotlib import pyplot as plt
from matplotlib import colors
from pyproj import Proj
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from string import replace
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator
from mudpy import gmttools
from mudpy.viewFQ import get_subfault_corners

faults_to_plot=5# 1 to 5main,synth,pian,anti,detach
#Run time parameters
home='/Users/dmelgar/Slip_inv/'

#project_name='Amatrice_M5.3-M5.8'
project_name='Amatrice_M6.6'
#project_name='Amatrice_complex'
#project_name='Amatrice_M6.1_August'
#project_name='Amatrice_M6.1_October's

run_name='gps_insar_vettore_main_synth_pian_anti_detach'
run_number='0022'

meshfile1='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_main.mshout'
meshfile2='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_synth.mshout'
meshfile3='/Users/dmelgar/Amatrice2016/3D_fault/_complex/pian_grande.mshout'
meshfile4='/Users/dmelgar/Amatrice2016/3D_fault/_complex/antithetic_listric.mshout'

#meshfile1='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_main.mshout'
#meshfile2='/Users/dmelgar/Amatrice2016/3D_fault/_complex/pian_grande.mshout'
#meshfile3='/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_synth.mshout'
#meshfile4='/Users/dmelgar/Amatrice2016/3D_fault/_complex/antithetic_listric.mshout'

af=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records.txt',usecols=[0,1,2,3,4,5])
af=array([[0,0,0,0,0]])
maxslip=None
UTM_zone='33T'
fudge=0.01
fake_hypo=[0.1,0.1],
borderwidth=0.5
figsize=(18,8)
xtick=0.2
ytick=0.2
ztick=2
inverse_model=False
hypocenter=None
strike=155
xl=[12.9,13.4]
yl=[42.5,43.2]
zl=[0,10]

whitejet_dict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.2, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
                'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.2, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
                'blue': ((0., 1, 1),
                  (0.05, 1, 1),
                  (0.2, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}


fault_name=home+project_name+'/output/inverse_models/models/%s.%s.inv' % (run_name,run_number)
gmttools.make_total_model(fault_name,thresh=0)
fault=genfromtxt(home+project_name+'/output/inverse_models/models/%s.%s.inv.total' % (run_name,run_number))
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
total_slip=slip

#Saturate to maxslip
if maxslip!=None:
    imax=where(slip>maxslip)[0]
    slip[imax]=maxslip
#normalize
norm_slip=slip/slip.max()

#get detachment corners
if faults_to_plot==5:
    nrows=4
    fault_detach=fault[len(corners1)+len(corners2)+len(corners3)+len(corners4):-1]
    fault_detach=fault_detach[0:-nrows*17+1]
    fault_detach=fault_detach[2*17:,:]
    slip_detach=norm_slip[len(corners1)+len(corners2)+len(corners3)+len(corners4):-1]
    corners_detach=get_subfault_corners(fault_detach)

#Get colormaps
#slip_colormap = colors.LinearSegmentedColormap('slip',whitejet_dict,256)
slip_colormap = plt.cm.gist_heat_r



#Plot init, axes positions etc
def make_subplot(ax1,sub_corners,sub_slip,af,azimuth,elevation,title,colorbar,triangles=True):
    
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
    ax1.scatter(af[:,3],af[:,2],af[:,4],s=5,lw=0)
    
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
    ax1.set_xlabel('\n\nLongitude')
    ax1.set_ylabel('\n\nLatitude')
    ax1.set_zlabel('Depth (km)',rotation=90)
    ax1.set_title(title)
    

#Determine whch indices for which fault
i_vettore_main=arange(0,len(corners1))
i_vettore_synth=arange(len(corners1),len(corners1)+len(corners2))
i_pian=arange(len(corners1)+len(corners2),len(corners1)+len(corners2)+len(corners3))
i_anti=arange(len(corners1)+len(corners2)+len(corners3),len(corners1)+len(corners2)+len(corners3)+len(corners4))

fig=plt.figure(figsize=figsize)

#Plot all faults
ax = fig.add_subplot(231, projection='3d')
if faults_to_plot==5:
    make_subplot(ax,corners_detach,slip_detach,af,azimuth=-60,elevation=51,title='',colorbar=False,triangles=False)
make_subplot(ax,corners,norm_slip,af,azimuth=-60,elevation=51,title='All',colorbar=False)

#Focus on Vettore main
ax = fig.add_subplot(232, projection='3d')
if faults_to_plot==1:
    colorbar=True
else:
    colorbar=False
make_subplot(ax,corners[i_vettore_main],norm_slip[i_vettore_main],af,azimuth=-100,elevation=57,title='Mt. Vettore main',colorbar=colorbar)

#Focus on vettore synth
if faults_to_plot>=2:
    ax = fig.add_subplot(233, projection='3d')
    if faults_to_plot==2  or faults_to_plot==5:
        colorbar=True
    else:
        colorbar=False
    make_subplot(ax,corners[i_vettore_synth],norm_slip[i_vettore_synth],af,azimuth=-100,elevation=57,title='Mt. Vettore synthetic',colorbar=colorbar)
    #make_subplot(ax,corners[i_vettore_synth],norm_slip[i_vettore_synth],af,azimuth=-10,elevation=46,title='Pian Grande',colorbar=colorbar)

#Focus on pian grande
if faults_to_plot>=3:
    if faults_to_plot==3:
        colorbar=True
    else:
        colorbar=False
    ax = fig.add_subplot(234, projection='3d')
    make_subplot(ax,corners[i_pian],norm_slip[i_pian],af,azimuth=-10,elevation=46,title='Pian Grande',colorbar=colorbar)

#Focus on antithetic
if faults_to_plot>=4:
    if faults_to_plot==4:
        colorbar=True
    else:
        colorbar=False
    ax = fig.add_subplot(235, projection='3d')
    make_subplot(ax,corners[i_anti],norm_slip[i_anti],af,azimuth=-10,elevation=46,title='Antithetic',colorbar=colorbar)

#Focus on detachment
if faults_to_plot==5:
    ax = fig.add_subplot(236, projection='3d')
    make_subplot(ax,corners_detach,slip_detach,af,azimuth=-60,elevation=51,title='Detachment',colorbar=True,triangles=False)

plt.subplots_adjust(hspace=0.15,wspace=0.1,left=0,right=0.95,top=0.99,bottom=0.07)

plt.show()