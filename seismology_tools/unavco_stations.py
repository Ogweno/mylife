from numpy import savetxt,array,c_
from matplotlib import pyplot as plt

f=open('/Users/dmelgar/Misc/UNAVCO2015.txt','r')

lonout=[]
latout=[]
while True:
    line=f.readline()
    if line=='':
        break
    latout.append(float(line.split('","')[3]))
    lonout.append(float(line.split('","')[4]))
lonout=array(lonout)
latout=array(latout)
plt.figure()
plt.scatter(lonout,latout)
savetxt('/Users/dmelgar/Misc/UNAVCO_parsed_2015.txt',c_[lonout,latout],fmt='%10.4f\t%10.4f')
plt.show()