from pyproj import Geod
from numpy import array,r_,ones,where
from string import rjust

f=open('/Users/dmelgar/Michoacan1985/gps/mexgps.txt','r')
fout=open('/Users/dmelgar/Michoacan1985/gps/gps.sta','w')
lonepi=-102.57
latepi=18.18
maxdist=400e3

k=0
while True:
    line=f.readline()
    if '>' in line:
        pass
    else:
        coords=array(line.split()).astype('float')
        if k==0:
            out=coords.copy()
        else:
            out=r_[out,coords]
        k+=1
    if line=='':
        break
f.close()
out=out.reshape((len(out)/2,2)).squeeze()   

#Select based on distance
p=Geod(ellps='WGS84')
az,baz,dist=p.inv(ones(len(out))*lonepi,ones(len(out))*latepi,out[:,0],out[:,1])
i=where(dist<maxdist)[0]
out=out[i,:]

for k in range(len(out)):
    sta='M'+rjust(str(k),3,'0')
    line='%s\t%.4f\t%.4f\n' % (sta,out[k,0],out[k,1])
    fout.write(line)
fout.close()
    
        
        
