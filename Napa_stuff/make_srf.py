from mudpy import forward
from glob import glob
#import bbp_plot
import os
from matplotlib import pyplot as plt

#logs=glob(u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/*.log')
#rupts=glob(u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/*.rupt')
#logs=[u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/napa.1sub.log']
#rupts=[u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/napa.1sub.rupt']
logs=glob(u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/napa.000000.log')
rupts=glob(u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/napa.000000.rupt')

for k in range(len(logs)):
    print k
    r=rupts[k]
    l=logs[k]
    #forward.mudpy2srf(r,l,time_pad=0,stf_type='dreger',stf_dt=0.01,minSTFpoints=24,integrate=True)
    forward.mudpy2sw4source(r,time_offset=0.0)
    
#srf=glob(u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/*.srf')
#src=glob(u'/Users/dmelgar/FakeQuakes/napa/output/ruptures/*.src')
##
#for k in range(len(srf)):
#    print k
#    fout='/Users/dmelgar/Napa2014/stochastic/plots/'+os.path.split(srf[k])[1].replace('.srf','.png')
#    bbp_plot.plot_srf(srf[k],src[k],1.0,xlim=[-15,6],figsize=[8,6],highlight=30,save=True,fout=fout)
#    
##plt.close('all')