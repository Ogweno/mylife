from shlex import split
import subprocess
from shutil import move,rmtree

base_path='/Users/dmelgar/code/BBP/bbp/bbp_data/'
#base_mw=['4.5','5.0','5.5','6.0','6.5','7.0']
base_mw=['5.5']
out_path='/Users/dmelgar/code/BBP/bbp/bbp_data/finished/large_eew_test/'
xml_path='/Users/dmelgar/code/BBP/bbp/bbp_data/run/xml/'
hypo_fractions=[0.05,0.25,0.5,0.75,0.95]

for kmag in range(len(base_mw)):
    for kfrac in range(len(hypo_fractions)):
        if kmag>10:
            pass
        #elif kmag==4 and kfrac>=1 and kfrac<3 :
        #    pass
        else:
            print 'Running m%s_frac%s' % (base_mw[kmag],hypo_fractions[kfrac])
            
            #Delete stuff first
            rmtree(base_path+'indata/123',ignore_errors=True)
            rmtree(base_path+'tmpdata/123',ignore_errors=True)  
            rmtree(base_path+'outdata/123',ignore_errors=True)
            rmtree(out_path+'m%s_frac%s' % (base_mw[kmag],hypo_fractions[kfrac]),ignore_errors=True)
            
            #Form xml file name for run
            xml_file=xml_path+'large_eew_%s_frac%s.xml' % (base_mw[kmag],hypo_fractions[kfrac])
            
            command='run_bbp.py -x %s -s 123' % (xml_file)
            print command
            command=split(command)
            p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err=p.communicate() 
            
            #Now move things/delete things
            move(base_path+'outdata/123/',out_path+'m%s_frac%s' % (base_mw[kmag],hypo_fractions[kfrac]))
        