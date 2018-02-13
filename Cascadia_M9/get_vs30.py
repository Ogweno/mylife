from numpy import genfromtxt,argmin

sta=genfromtxt(u'/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl')
sta_name=genfromtxt(u'/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl',usecols=[2],dtype='S')
vs30=genfromtxt(u'/Users/dmelgar/DEMs/vs30/cascadia_vs30.xyz')

for k in range(len(sta)):
    dist=((sta[k,0]-vs30[:,0])**2+(sta[k,1]-vs30[:,1])**2)**0.5
    i=argmin(dist)
    vs30_sta=vs30[i,2]
    print sta_name[k]+'\t'+str(int(vs30_sta))