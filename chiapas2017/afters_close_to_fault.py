from matplotlib import pyplot as plt
from numpy import genfromtxt,where,isnan,argmax,ones,array,r_,c_,savetxt,expand_dims
from pyproj import Geod

afters_dist=30*1e3
afters=genfromtxt('/Users/dmelgar/Chiapas2017/afters/SSNMX_cleaned.csv',delimiter=',',usecols=[3,4,5])
fault=genfromtxt('/Users/dmelgar/Chiapas2017/slab/chiapas_normal_ssn_trim.fault')

p=Geod(ellps='WGS84')
#Get aftershocks within some distance of fault
af=array([])
count=0
for k in range(len(afters)):
    az,baz,dist=p.inv(fault[:,1],fault[:,2],afters[k,1]*ones(len(fault)),afters[k,0]*ones(len(fault)))
    if dist.min()<afters_dist:
        if count==0:
            af=expand_dims(afters[k,:],0)
            count+=1
        else:    
            af=r_[af,expand_dims(afters[k,:],0)]
af[:,2]=-af[:,2]
af_out=af.copy()
af_out[:,0]=af[:,1]
af_out[:,1]=af[:,0]
af_out[:,2]=af[:,2]
savetxt('/Users/dmelgar/Chiapas2017/afters/afters_10days_close2fault.txt',af_out,fmt='%.4f')