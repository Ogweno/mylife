from numpy import genfromtxt,where,float64,r_,expand_dims
from glob import glob
from gmsh_tools import xyz2gmsh

#Gmsh output file name
gmsh_out=u'/Users/dmelgar/Melinka2016/etc/slab3D/melinka.gmsh'
ztop=-19
zbot=-45


#Line filters
[L1x1 , L1y1] = [-74.6065 , -43.8669]
[L1x2 , L1y2] = [-73.6540 , -44.0468]
[L2x1 , L2y1] = [-74.5505 , -42.8613]
[L2x2 , L2y2] = [-73.4977 , -42.7315]

#Contours files 
contour_files= glob('/Users/dmelgar/Melinka2016/etc/slab3D/*.contour')
#read fault
f=genfromtxt('/Users/dmelgar/Melinka2016/etc/slab_tasa.xyz')
#Read top limit
perimeter=genfromtxt(u'/Users/dmelgar/Melinka2016/perimeter.xy')

#Parse contours file
#contours=top.copy()
first=True
for k in range(len(contour_files)):
    f=open(contour_files[k])
    while True:
        line=f.readline().rstrip()
        if not line: #Exit loop
            break
        #Assign line info to array
        if '>' not in line: #It's not a segment header
            if first==True:
                contours=expand_dims(float64(line.split('\t')),0)
                first=False
            else:
                contours=r_[contours,expand_dims(float64(line.split('\t')),0)]   
    f.close()

#Now filter things
#By depth
i=where(contours[:,2]>=zbot)[0]
contours=contours[i,:]
i=where(contours[:,2]<=ztop)[0]
contours=contours[i,:]
#By regions
#Equations of the line filters
m1=(L1y1-L1y2)/(L1x1-L1x2)
m2=(L2y1-L2y2)/(L2x1-L2x2)
b1=L1y1-m1*L1x1
b2=L2y1-m2*L2x1
#Get points ABOVE line 1
ytest=m1*contours[:,0]+b1
i=where((ytest-contours[:,1])<=0)[0] #Points that lie below the line
contours=contours[i,:]
#Get points BELOW line 2
ytest=m2*contours[:,0]+b2
i=where((ytest-contours[:,1])>=0)[0] #Points that lie below the line
contours=contours[i,:]

#Write gmsh file
xyz2gmsh(gmsh_out,contours[:,0],contours[:,1],contours[:,2],coord_type='UTM')


    
    




