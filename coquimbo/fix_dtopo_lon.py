from glob import glob
from numpy import genfromtxt,savetxt

#folders=glob('/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule*')
folders=['/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_7.6930.sub0022',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.3260.sub0032',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.4310.sub0175',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.4510.sub0048',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.4960.sub0044',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.5280.sub0155',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.5690.sub0063',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.5690.sub0090',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.5940.sub0088',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.6090.sub0079',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.6500.sub0096',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_8.7120.sub0071',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.0050.sub0029',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.2690.sub0178',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.4020.sub0085',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.4100.sub0030',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.4130.sub0186',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.4350.sub0098',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.5090.sub0042',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.5430.sub0180',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.5820.sub0102',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.6350.sub0093',
'/Users/dmelgar/Slip_inv/Coquimbo/GFs/tsunami/maule.mod_9.6700.sub0129']

for k in range(len(folders)):
    print folders[k]
    ds=genfromtxt(folders[k]+'/DS.dtopo')
    ss=genfromtxt(folders[k]+'/SS.dtopo')
    if ds[1,1]>0:
        ds[:,1]=ds[:,1]-360
        ss[:,1]=ss[:,1]-360
    savetxt(folders[k]+'/DS.dtopo',ds,fmt='%i\t%12.6f\t%12.6f\t%.6e')
    savetxt(folders[k]+'/SS.dtopo',ss,fmt='%i\t%12.6f\t%12.6f\t%.6e')