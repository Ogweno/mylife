from glob import glob
from os import remove

dirs=glob('/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/tsunami/maule*')

for k in range(len(dirs)):
    print dirs[k]
    files=glob(dirs[k]+'/'+'*grn*')
    for j in range(len(files)):
        remove(files[j])
    