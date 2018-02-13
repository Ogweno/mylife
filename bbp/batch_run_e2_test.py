from shlex import split
import subprocess
from shutil import move,rmtree
from numpy.random import rand

base_path='/Users/dmelgar/code/BBP/bbp/bbp_data/'
base_mw='5.0'
iterations=50
out_path='/Users/dmelgar/code/BBP/bbp/bbp_data/finished/e2_test/'
xml_path='/Users/dmelgar/code/BBP/bbp/bbp_data/run/xml/'
xml_file=xml_path+'e2_m5.0.xml'
src_path='/Users/dmelgar/code/BBP/bbp/bbp_data/run/src/e2_test/'
src_template='m5.0_template.src'
src_out='m5.0.src'

for k in range(iterations):

    print 'Running iteration %d' % (k)
    
    #Delete stuff first
    rmtree(base_path+'indata/123',ignore_errors=True)
    rmtree(base_path+'tmpdata/123',ignore_errors=True)  
    rmtree(base_path+'outdata/123',ignore_errors=True)
    rmtree(out_path+'m%s_iteration%d' % (base_mw,k),ignore_errors=True)
    
    #Change SEED value
    f=open(src_path+src_template,'r')
    fwrite=open(src_path+src_out,'w')
    while True:
        line=f.readline()
        if line=='':
            break
        if 'SEED' in line:
            seed=int(rand()*1e6)
            line='SEED = '+str(seed)
        fwrite.write(line)
    f.close()
    fwrite.close()
    
    
    command='run_bbp.py -x %s -s 123' % (xml_file)
    print command
    command=split(command)
    p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err=p.communicate() 
    
    #Now move things/delete things
    move(base_path+'outdata/123/',out_path+'m%s_iteration%d' % (base_mw,k))
