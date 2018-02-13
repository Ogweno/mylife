from numpy import genfromtxt,r_,zeros,unique,where,diff
from matplotlib import pyplot as plt

variable='Mzz'
fout=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished1.los'
f0=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process0'
f1=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process1'
f2=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process2'
f3=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process3'
f4=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process4'
f5=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process5'
f6=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process6'
f7=u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process7'

filenames=[f0,f1,f2,f3,f4,f5,f6,f7]

with open(fout, 'w') as outfile:
    outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
    for fname in filenames:
        print fname
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
            outfile.write('\n')
                