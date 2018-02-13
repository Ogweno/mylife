from numpy import genfromtxt,c_,r_,savetxt,array
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

msh=genfromtxt('/Users/dmelgar/Melinka2016/etc/slab3D/melinka.mshout')
f=genfromtxt('/Users/dmelgar/Slip_inv/Melinka_3D/output/inverse_models/models/jacknife.0174.inv')
f=f[0:365,:]
fout='/Users/dmelgar/Melinka2016/etc/slab3D/checker.rupt'
select=[61,55,35,14,8,1,5,11,30,46,19,16,6,66,47,18,43,31,34,3,29,4,62,2,36,21,0,25,32,
    122,128,119,143,138,127,75,91,145,136,129,123,73,85,113,98,71,101,113,112,70,74,150,
    144,111,103,84,106,115,69,78,87,132,99,96,125,118,88,89,82,117,124,109,95,97,83,108,133,148,135,134,
    159,155,228,153,154,222,202,199,206,201,200,185,
    160,204,183,229,174,190,192,225,158,172,205,189,184,214,173,188,212,197,170,
    287,243,251,252,265,256,258,279,270,257,268,274,268,271,235,247,286,278,259,242,255,
    263,262,276,241,275,244,236,289,234,266,253,285,269,282,277,288,245,250,
    294,347,350,358,363,321,317,315,329,307,334,311,343,333,
    340,345,344,359,301,304,308,303,361,341,323,353,326,336,353,342,339,337]

plt.close('all')
fig=plt.figure(figsize=(16,16))
ax=fig.add_subplot(111)
patches=[]
for k in range(len(msh)):
    node1=msh[k,4:7]
    node2=msh[k,7:10]
    node3=msh[k,10:13]
    c=msh[k,1:3]
    c[0]=c[0]+360
    c[0]=c[0]-0.02
    l=str(k)
    x=r_[node1[0],node2[0],node3[0],node1[0]]
    y=r_[node1[1],node2[1],node3[1],node1[1]]
    verts=c_[x,y]
    plt.plot(x,y,'k')
    plt.annotate(l,xy=c)
    if k in select:
        polygon = Polygon(verts,closed=True)
        patches.append(polygon)

collection = PatchCollection(patches)
ax.add_collection(collection)
collection.set_color('r')

plt.axis('equal')        
plt.show()

#write file
#      1	 -74.4006	 -43.1867	 20.4800	 352.24	15.04	0.50	2.00	  2.0636e-01	  2.0636e-01	    6708.9	    6708.9	 17.6335	 4.4790e+10
f[:,7]=4.0
f[:,8]=0
f[:,9]=0
select=array(select)
f[select,9]=1.0
savetxt(fout,f,fmt='%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%.1f\t%.1f\t%.4f\t%.4e')