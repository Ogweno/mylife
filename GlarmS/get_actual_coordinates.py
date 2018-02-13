from numpy import genfromtxt,where,zeros

f=open(u'/Users/dmelgar/GlarmS/GlarmS_west_coast_stations.txt','w')
names=genfromtxt(u'/Users/dmelgar/GlarmS/actaully_used_list.txt',dtype='S')
all_names=genfromtxt(u'/Users/dmelgar/GlarmS/sites_RT_02082017.list',usecols=3,dtype='S')
all_cords=genfromtxt(u'/Users/dmelgar/GlarmS/sites_RT_02082017.list',usecols=[0,1])
coords=zeros((len(names),2))
for k in range(len(names)):
    i=where(all_names==names[k])[0][0]
    print names[k],i
    coords[k,:]=all_cords[i,:]
    
#now save
for k in range(len(coords)):
    line='%s\t%.6f\t%.6f\n' % (names[k],coords[k,0],coords[k,1])
    f.write(line)
f.close()