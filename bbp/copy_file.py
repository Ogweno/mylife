from shutil import copy
from string import rjust
path=u'/Users/dmelgar/FakeQuakes/M6_validation/output/ruptures/'

for k in range(1,100):
    name='M6.'+rjust(str(k),6,'0')+'.log'
    copy( u'/Users/dmelgar/FakeQuakes/M6_validation/output/ruptures/M6.000000.log',path+name)
    
    name='M6.'+rjust(str(k),6,'0')+'.rupt'
    copy( u'/Users/dmelgar/FakeQuakes/M6_validation/output/ruptures/M6.000000.rupt',path+name)