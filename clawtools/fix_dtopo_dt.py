from glob import glob
from numpy import genfromtxt,savetxt

dirs=glob(u'/Users/dmelgar/Slip_inv/iquique_joint/GFs/tsunami/*sub*')
correct=4

for k in range(len(dirs)):
    print k
    DS=genfromtxt(dirs[k]+'/DS.dtopo')
    SS=genfromtxt(dirs[k]+'/SS.dtopo')
    DS[:,0]=DS[:,0]*4
    SS[:,0]=SS[:,0]*4
    savetxt(dirs[k]+'/DS.dtopo',DS,fmt='%d\t%.6f\t%.6f\t%.6e')
    savetxt(dirs[k]+'/SS.dtopo',SS,fmt='%d\t%.6f\t%.6f\t%.6e')