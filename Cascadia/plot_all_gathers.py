from mudpy import viewFQ
from glob import glob
from matplotlib import pyplot as plt

component='n'
home='/Users/dmelgar/fakequakes/'
project_name='Cascadia'
GF_list='cascadia_small.gflist'
out_folder='/Users/dmelgar/FakeQuakes/Cascadia/plots/_gathers/'
ruptures=glob('/Users/dmelgar/FakeQuakes/Cascadia/output/ruptures/*.rupt')

for k in range(len(ruptures)):
    print k
    current_rupture=ruptures[k].split('/')[-1].split('.')[0]+'.'+ruptures[k].split('/')[-1].split('.')[1]
    viewFQ.record_section(home,project_name,GF_list,current_rupture,factor=40)
    plt.savefig(out_folder+current_rupture+'.gather.png')
    plt.close()

