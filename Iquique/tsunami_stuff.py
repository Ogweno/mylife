'''
Iquique tsunamis tuff
'''
from numpy import arange,zeros,savetxt,c_,where
from string import rjust

make_grid=True
gridout='/Users/dmelgar/Maule2010/tsunami/tsun_grid.xy'

lon=[-76.2,-71]
lat=[-38.85,-33]
delta=0.075

m1=(-38.4946+30.9728)/(-76.3853+73.5890)
b1=-38.4946-m1*(-76.3853)
m2=(-39.4400+30.8607)/(-73.0604+70.9900)
b2=-39.44-m2*(-73.0604)
m3=(-37.7414+38.9361)/(-76.7110+71.9067)
b3=-37.7414-m3*(-76.7110)
m4=(-32.0578+33.152)/(-74.4251+70.3117)
b4=-32.0578-m4*(-74.4251)
if make_grid==True:
    longrid=arange(lon[0],lon[1]+0.0001,delta)
    latgrid=arange(lat[0],lat[1]+0.0001,delta)
    latout=zeros(len(longrid)*len(latgrid))
    lonout=zeros(len(longrid)*len(latgrid))
    k=0
    for i in range(len(latgrid)):
        for j in range(len(longrid)):
            latout[k]=latgrid[-(i+1)]
            lonout[k]=longrid[j]
            k+=1
    #Now filter based on lines
    ytest=m1*lonout+b1
    i=where(ytest>latout)[0]
    lonout=lonout[i]
    latout=latout[i]
    ytest=m2*lonout+b2
    i=where(ytest<latout)[0]
    lonout=lonout[i]
    latout=latout[i]
    ytest=m3*lonout+b3
    i=where(ytest<latout)[0]
    lonout=lonout[i]
    latout=latout[i]
    ytest=m4*lonout+b4
    i=where(ytest>latout)[0]
    lonout=lonout[i]
    latout=latout[i]
    f=open(gridout,'w')
    for k in range(len(lonout)):
        line='tg%s\t%10.6f\t%10.6f\n' %(rjust(str(k),4,'0'),lonout[k],latout[k])
        f.write(line)
    f.close()
    