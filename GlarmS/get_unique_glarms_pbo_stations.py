from numpy import genfromtxt,r_,array,expand_dims

stations_glarms=genfromtxt('/Users/dmelgar/GlarmS/GlarmS_west_coast_stations.txt',usecols=[1,2])
names_glarms=genfromtxt('/Users/dmelgar/GlarmS/GlarmS_west_coast_stations.txt',usecols=0,dtype='S')
stations_pbo=genfromtxt(u'/Users/dmelgar/GlarmS/list_newpbo.txt',usecols=[2,1])
names_pbo=genfromtxt(u'/Users/dmelgar/GlarmS/list_newpbo.txt',usecols=[0],dtype='S')

stations_out=stations_pbo.copy()
names_out=names_pbo.copy()
fout='/Users/dmelgar/GlarmS/Glarm_pbo_merge.txt'
#Loop over glarms stations if no close station then ADD to pbo list
for k in range(len(stations_glarms)):
    
    dist=((stations_glarms[k,0]-stations_pbo[:,0])**2+(stations_glarms[k,1]-stations_pbo[:,1])**2)**0.5
    
    if dist.min()<0.001: #Duplicate DON'T ADD
        print '... IGNORING: GlarmS name is '+str(names_glarms[k])+' ; PBO name is '+str(names_pbo[dist.argmin()])+' ; dist is'+str(dist.min())
        
    else:  #NOT a duplicate, ADD
        print '... ADDING: GlarmS name is '+str(names_glarms[k])
        names_out=r_[names_out,array([names_glarms[k]])]
        stations_out=r_[stations_out,expand_dims(stations_glarms[k,:],0)]

f=open(fout,'w')
for k in range(len(names_out)):
    line='%s\t%.6f\t%.6f\n' % (names_out[k],stations_out[k,0],stations_out[k,1])
    f.write(line)
f.close()