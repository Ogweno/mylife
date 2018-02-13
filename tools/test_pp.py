import pp
from numpy import loadtxt
import numpy
import time

def printsource(s,k):
    #from time import sleep
    #from numpy.random import rand
    t=numpy.random.rand()*1
    time.sleep(t)
    print 'Source number',k
    print s
    print 'Slept for',t,'seconds'

ncpus=1
fault_file=u'/Users/dmelgar/Slip_inv/paralell/data/model_info/nepal_10.fault'
source=loadtxt(fault_file,ndmin=2)
ppservers = ()
jobs=[]
job_server = pp.Server(ncpus, ppservers=ppservers)
print "Starting GFs computation on", job_server.get_ncpus(), "CPUs"
for k in range(len(source)):
    jobs.append(job_server.submit(printsource,(source[k,:],k),(),("numpy","time")))
    
#for job in jobs:
#    result = job()
    
job_server.print_stats()
job_server.destroy()