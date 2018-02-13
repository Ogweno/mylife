from pyproj import Geod
from numpy import array,zeros,sqrt,arange
from numpy.random import rand
from string import rjust

#f=open('/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/eew_test.stl','w')
#f=open('/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/fake_nocal.stl','w')
#f=open('/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/pwave.stl','w')
f=open('/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/e2.stl','w')
#f2=open('/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/fake_nocal.channels.txt','w')

hypo=array([-122.26,37.87])
g=Geod(ellps='WGS84')

#lon=zeros(46)
#lat=zeros(46)



##Station coords
#lon[0],lat[0],baz=g.fwd(hypo[0],hypo[1],270,5000)
#lon[1],lat[1],baz=g.fwd(hypo[0],hypo[1],270,10000)
#lon[2],lat[2],baz=g.fwd(hypo[0],hypo[1],270,15000)
#
#lon[3],lat[3],baz=g.fwd(lon[0],lat[0],0,5000)
#lon[4],lat[4],baz=g.fwd(lon[0],lat[0],0,10000)
#lon[5],lat[5],baz=g.fwd(lon[0],lat[0],180,5000)
#lon[6],lat[6],baz=g.fwd(lon[0],lat[0],180,10000)
#
#lon[7],lat[7],baz=g.fwd(lon[3],lat[3],270,5000)
#lon[8],lat[8],baz=g.fwd(lon[3],lat[3],270,10000)
#
#lon[9],lat[9],baz=g.fwd(lon[4],lat[4],270,5000)
#lon[10],lat[10],baz=g.fwd(lon[4],lat[4],270,10000)
#
#lon[11],lat[11],baz=g.fwd(lon[5],lat[5],270,5000)
#lon[12],lat[12],baz=g.fwd(lon[5],lat[5],270,10000)
#
#lon[13],lat[13],baz=g.fwd(lon[6],lat[6],270,5000)
#lon[14],lat[14],baz=g.fwd(lon[6],lat[6],270,10000)



#P-waves
#lon[0],lat[0],baz=g.fwd(hypo[0],hypo[1],45,5000)
#lon[1],lat[1],baz=g.fwd(hypo[0],hypo[1],45,10000)
#lon[2],lat[2],baz=g.fwd(hypo[0],hypo[1],45,20000)
#lon[3],lat[3],baz=g.fwd(hypo[0],hypo[1],45,50000)
#lon[4],lat[4],baz=g.fwd(hypo[0],hypo[1],45,100000)
#lon[5],lat[5],baz=g.fwd(hypo[0],hypo[1],45,200000)
#lon[6],lat[6],baz=g.fwd(hypo[0],hypo[1],45,400000)

#Regular grid
longrid=arange(hypo[0]-0.5,hypo[0]+0.51,0.1)
latgrid=arange(hypo[1]-0.5,hypo[1]+0.51,0.1)
k=0
lon=zeros(len(longrid)*len(latgrid))
lat=zeros(len(longrid)*len(latgrid))
for i in range(len(longrid)):
    for j in range(len(latgrid)):
        lon[k]=longrid[i]
        lat[k]=latgrid[j]
        k+=1



#Ciruclar

##Station coords
#lon[0],lat[0],baz=g.fwd(hypo[0],hypo[1],270,sqrt(1000**2+1000**2))
#
#lon[1],lat[1],baz=g.fwd(lon[0],lat[0],0,1000)
#lon[2],lat[2],baz=g.fwd(lon[0],lat[0],0,5000)
#lon[3],lat[3],baz=g.fwd(lon[0],lat[0],0,10000)
#lon[4],lat[4],baz=g.fwd(lon[0],lat[0],0,20000)
#lon[5],lat[5],baz=g.fwd(lon[0],lat[0],0,40000)
#
#lon[6],lat[6],baz=g.fwd(lon[0],lat[0],270,1000)
#lon[7],lat[7],baz=g.fwd(lon[0],lat[0],270,5000)
#lon[8],lat[8],baz=g.fwd(lon[0],lat[0],270,10000)
#lon[9],lat[9],baz=g.fwd(lon[0],lat[0],270,20000)
#lon[10],lat[10],baz=g.fwd(lon[0],lat[0],270,40000)
#
#lon[11],lat[11],baz=g.fwd(lon[0],lat[0],292.5,1000)
#lon[12],lat[12],baz=g.fwd(lon[0],lat[0],292.5,5000)
#lon[13],lat[13],baz=g.fwd(lon[0],lat[0],292.5,10000)
#lon[14],lat[14],baz=g.fwd(lon[0],lat[0],292.5,20000)
#lon[15],lat[15],baz=g.fwd(lon[0],lat[0],292.5,40000)
#
#lon[16],lat[16],baz=g.fwd(lon[0],lat[0],315,1000)
#lon[17],lat[17],baz=g.fwd(lon[0],lat[0],315,5000)
#lon[18],lat[18],baz=g.fwd(lon[0],lat[0],315,10000)
#lon[19],lat[19],baz=g.fwd(lon[0],lat[0],315,20000)
#lon[20],lat[20],baz=g.fwd(lon[0],lat[0],315,40000)
#
#lon[21],lat[21],baz=g.fwd(lon[0],lat[0],337.5,1000)
#lon[22],lat[22],baz=g.fwd(lon[0],lat[0],337.5,5000)
#lon[23],lat[23],baz=g.fwd(lon[0],lat[0],337.5,10000)
#lon[24],lat[24],baz=g.fwd(lon[0],lat[0],337.5,20000)
#lon[25],lat[25],baz=g.fwd(lon[0],lat[0],337.5,40000)
#
#lon[26],lat[26],baz=g.fwd(lon[0],lat[0],247.5,1000)
#lon[27],lat[27],baz=g.fwd(lon[0],lat[0],247.5,5000)
#lon[28],lat[28],baz=g.fwd(lon[0],lat[0],247.5,10000)
#lon[29],lat[29],baz=g.fwd(lon[0],lat[0],247.5,20000)
#lon[30],lat[30],baz=g.fwd(lon[0],lat[0],247.5,40000)
#
#lon[31],lat[31],baz=g.fwd(lon[0],lat[0],225,1000)
#lon[32],lat[32],baz=g.fwd(lon[0],lat[0],225,5000)
#lon[33],lat[33],baz=g.fwd(lon[0],lat[0],225,10000)
#lon[34],lat[34],baz=g.fwd(lon[0],lat[0],225,20000)
#lon[35],lat[35],baz=g.fwd(lon[0],lat[0],225,40000)
#
#lon[36],lat[36],baz=g.fwd(lon[0],lat[0],202.5,1000)
#lon[37],lat[37],baz=g.fwd(lon[0],lat[0],202.5,5000)
#lon[38],lat[38],baz=g.fwd(lon[0],lat[0],202.5,10000)
#lon[39],lat[39],baz=g.fwd(lon[0],lat[0],202.5,20000)
#lon[40],lat[40],baz=g.fwd(lon[0],lat[0],202.5,40000)
#
#lon[41],lat[41],baz=g.fwd(lon[0],lat[0],180,1000)
#lon[42],lat[42],baz=g.fwd(lon[0],lat[0],180,5000)
#lon[43],lat[43],baz=g.fwd(lon[0],lat[0],180,10000)
#lon[44],lat[44],baz=g.fwd(lon[0],lat[0],180,20000)
#lon[45],lat[45],baz=g.fwd(lon[0],lat[0],180,40000)



#CRandom

##Station coords
#
#for k in range(len(lon)):
#    radius=rand(1)*20000
#    azimuth=rand(1)*360
#    lon[k],lat[k],baz=g.fwd(hypo[0],hypo[1],azimuth,radius)
#
#
##JW  TEST -- CNZ 37.8741 -122.2598    40.0   100.00   416.666                 DU/M/S**2
##columns: network station location channel latitude longitude elevation samplerate gain units
#
#for k in range(len(lon)):
#    line='%.4f\t%.4f\tst%s\t720\n' % (lon[k],lat[k],rjust(str(k),2,'0'))
#    line2='DM\tst%s\t--\tBNZ\t%.4f\t%.4f\t0.0\t40.00\t0.01\tDU/M/S**2\n' % (rjust(str(k),2,'0'),lat[k],lon[k])
#    f.write(line)
#    f2.write(line2)
#f.close()
#f2.close()

for k in range(len(lon)):
    line='%.4f\t%.4f\tS%s\t720\n' % (lon[k],lat[k],rjust(str(k),3,'0'))
    f.write(line)
f.close()