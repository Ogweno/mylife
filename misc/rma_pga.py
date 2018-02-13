from numpy import genfromtxt,ones,zeros,savetxt,c_
import gmpe_tools
from pyproj import Geod

lonlat=genfromtxt('/Users/dmelgar/Downloads/20160718-query.csv',usecols=[1,2],delimiter=',')
m=genfromtxt('/Users/dmelgar/Downloads/20160718-query.csv',usecols=[4],delimiter=',')

g=Geod(ellps='WGS84')
lonsta=-122.9090*ones(len(m))
latsta=39.1060*ones(len(m))

az,baz,dist=g.inv(lonsta,latsta,lonlat[:,1],lonlat[:,0])
dist=dist/1000

vs30=720*ones(len(m))

U=zeros(len(m))
NS=zeros(len(m))
RS=zeros(len(m))
SS=ones(len(m))

pga,stdpga=gmpe_tools.bssa14(m,dist,vs30,intensity_measure='PGA',U=U,NS=NS,RS=RS,SS=SS)
pgv,stdpgv=gmpe_tools.bssa14(m,dist,vs30,intensity_measure='PGV',U=U,NS=NS,RS=RS,SS=SS)

savetxt('/Users/dmelgar/Downloads/RMA_groundmotion.txt',c_[m,dist,pga,stdpga,pgv,stdpgv],fmt='%.2f\t%.2f\t%.4e\t%.4f\t%.4e\t%.4f',header='M, dist(km), PGA (g), sigma_PGA, PGV (cm/s), sigma_PGV')