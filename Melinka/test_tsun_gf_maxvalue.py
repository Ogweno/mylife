from obspy import read
from glob import glob
from os import chdir
import subprocess
from shlex import split

p="/Users/dmelgar/Slip_inv/Melinka_usgs/GFs/tsunami/"
folders=glob(p+'*maule*')
thresh=1.0
count=0
for k in range(len(folders)):
    ds=read(folders[k]+'/pmel.ds.tsun')
    ss=read(folders[k]+'/pmel.ss.tsun')
    maxds=abs(ds[0].data.max())
    maxss=abs(ss[0].data.max())
    if maxds>thresh or maxss>thresh:
        count+=1
        print folders[k]
        print '... ds max is %f' %(maxds)
        print '... ss max is %f' %(maxss)
        ##ok go into_DS and _SS and modify setrun and run
        #s_ds = open(folders[k]+'/_DS/setrun.py').read()
        #s_ss = open(folders[k]+'/_SS/setrun.py').read()
        #s_ds = s_ds.replace('refinement_data.wave_tolerance = 0.001','refinement_data.wave_tolerance = 0.0005')
        #s_ss = s_ss.replace('refinement_data.wave_tolerance = 0.001','refinement_data.wave_tolerance = 0.0005')
        #fpy = open(folders[k]+'/_DS/setrun.py', 'w')
        #fpy.write(s_ds)
        #fpy.close()
        #fpy = open(folders[k]+'/_SS/setrun.py', 'w')
        #fpy.write(s_ss)
        #fpy.close()
        ##no run it
        #geoclawDS='make .output'
        #geoclawSS='make .output'
        #chdir(folders[k]+'/_DS')
        #subprocess.call(split('make clean'))
        #subprocess.call(split(geoclawDS))
        #chdir(folders[k]+'/_SS')
        #subprocess.call(split('make clean'))
        #subprocess.call(split(geoclawSS))
        
        
print '%d subs over threshold' %(count)