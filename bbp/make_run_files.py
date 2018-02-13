import bbptools
from numpy.random import rand


#base_mw=['4.5','5.0','5.5','6.0','6.5','7.0']
#xml_mags=['45','50','55','60','65','70']
#base_mw=['4.5','5.0','5.5','6.0','6.5']
#xml_mags=['45','50','55','60','65']
base_mw=['7.0']
xml_mags=['70']
base_filename='/Users/dmelgar/code/BBP/bbp/bbp_data/run/src/large_eew/template_m'
base_out_file='/Users/dmelgar/code/BBP/bbp/bbp_data/run/src/large_eew/large_eew_m'
base_xml='/Users/dmelgar/code/BBP/bbp/bbp_data/run/xml/'
hypo_fractions=[0.05,0.25,0.5,0.75,0.95]

make_src=False
make_xml=True

for kmag in range(len(base_mw)):
    for kfrac in range(len(hypo_fractions)):
        
        if make_src:
            base_src_file=base_filename+base_mw[kmag]+'.src'
            out_src_file=base_out_file+base_mw[kmag]+'_frac'+str(hypo_fractions[kfrac])+'.src'
            hypo_position=bbptools.hypo_position(base_src_file,hypo_fraction=hypo_fractions[kfrac])
            
            #Now make new file
            fin=open(base_src_file,'r')
            fout=open(out_src_file,'w')
            while True:
                line=fin.readline()
                if 'HYPO_ALONG_STK' in line:
                    fout.write('HYPO_ALONG_STK = %.4f\n' % (hypo_position))
                elif 'SEED' in line:
                    seed=int(rand(1)*100000)
                    fout.write('SEED = %d\n' % (seed))
                else:
                    fout.write(line)
                if line=='':
                    break
            fin.close()
            fout.close()
        
        #Now deal with xml
        if make_xml:
            xml_in=base_xml+xml_mags[kmag]+'.xml'
            xml_out=base_xml+'large_eew_'+base_mw[kmag]+'_frac'+str(hypo_fractions[kfrac])+'.xml'
            fin=open(xml_in,'r')
            fout=open(xml_out,'w')
            while True:
                line=fin.readline()
                test_src='template_m'+base_mw[kmag]+'.src'
                test_srf='template_m'+base_mw[kmag]+'.srf'
                test_path='/test_eew/'
                if test_src in line:
                    lineout=line.replace(test_src,'large_eew_m'+base_mw[kmag]+'_frac'+str(hypo_fractions[kfrac])+'.src')
                    if test_path in lineout:
                        lineout=lineout.replace(test_path,'/large_eew/')
                    fout.write(lineout)
                elif test_srf in line:
                    lineout=line.replace(test_srf,'large_eew_m'+base_mw[kmag]+'_frac'+str(hypo_fractions[kfrac])+'.srf')
                    fout.write(lineout)
                else:
                    fout.write(line)
                if line=='':
                    break
            fin.close()
            fout.close()
            
        