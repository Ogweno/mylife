from numpy import r_,array,c_,savetxt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


fseis='/Users/dmelgar/Coquimbo2015/afters/select.out'
fout='/Users/dmelgar/Coquimbo2015/afters/csn_seismicity.txt'

f=open(fseis,'r')
previous_line=''
lon=array([])
lat=array([])
depth=array([])
mag=array([])
dayout=array([])
monthout=array([])
yearout=array([])
while True:
    line=f.readline()
    try:
        if (line[-2] == '1'):
        #if 'GAP=' in line:
            day=float(line[8:10])
            mo=float(line[5:8])
            yr=float(line[0:5])
            if yr==2015:
                if mo<9:
                    print float(line[56:59])
                    lon=r_[lon,float(line[31:38])]
                    lat=r_[lat,float(line[23:30])]
                    depth=r_[depth,float(line[38:43])]
                    mag=r_[mag,float(line[56:59])]
                    dayout=r_[dayout,day]
                    monthout=r_[monthout,mo]
                    yearout=r_[yearout,yr]
            else:
                    print float(line[56:59])
                    lon=r_[lon,float(line[31:38])]
                    lat=r_[lat,float(line[23:30])]
                    depth=r_[depth,float(line[38:43])]
                    mag=r_[mag,float(line[56:59])]
                    dayout=r_[dayout,day]
                    monthout=r_[monthout,mo]
                    yearout=r_[yearout,yr]
    except:
        print line
    previous_line=line
    if line=='':
        break
out=c_[yearout,monthout,dayout,lon,lat,depth,mag]
savetxt(fout,out,fmt='%i\t%i\t%i\t%.4f\t%.4f\t%6.1f\t%6.1f',header='year,month,day,lon,lat,depth(km),magnitude')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(lon,lat,-depth)
plt.show()

        