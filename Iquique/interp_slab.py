from numpy import genfromtxt,where,float64,r_,expand_dims
from glob import glob
from gmsh_tools import xyz2gmsh

#Gmsh output file name
gmsh_out=u'/Users/dmelgar/Iquique2014/iquique.gmsh'

#Depth filter
maxdepth=-70
#Line filters
[L1x1 , L1y1] = [360-73.2593 , -19.3565]
[L1x2 , L1y2] = [360-68.5669 , -17.5311]
[L2x1 , L2y1] = [360-72.9443 , -21.0565]
[L2x2 , L2y2] = [360-68.4823 , -21.0778]

#Contours files 
contour_files= glob('/Users/dmelgar/Models/*.contour')
#read fault
f=genfromtxt('/Users/dmelgar/Models/sam_slab1.0_clip.xyz')
#Read top limit
top=genfromtxt('/Users/dmelgar/Models/sam_top.in')

#Parse contours file
contours=top.copy()
for k in range(len(contour_files)):
    f=open(contour_files[k])
    while True:
        line=f.readline().rstrip()
        if not line: #Exit loop
            break
        #Assign line info to array
        if '>' not in line: #It's not a segment header
            contours=r_[contours,expand_dims(float64(line.split('\t')),0)]   
    f.close()

#Now filter things
#By depth
i=where(contours[:,2]>=maxdepth)[0]
contours=contours[i,:]
#By regions
#Equations of the line filters
m1=(L1y1-L1y2)/(L1x1-L1x2)
m2=(L2y1-L2y2)/(L2x1-L2x2)
b1=L1y1-m1*L1x1
b2=L2y1-m2*L2x1
#Get points below line 1
ytest=m1*contours[:,0]+b1
i=where((ytest-contours[:,1])>=0)[0] #Points that lie below the line
contours=contours[i,:]
#Get points above line 2
ytest=m2*contours[:,0]+b2
i=where((ytest-contours[:,1])<=0)[0] #Points that lie below the line
contours=contours[i,:]

#Write gmsh file
xyz2gmsh(gmsh_out,contours[:,0],contours[:,1],contours[:,2],coord_type='UTM')


    
    




