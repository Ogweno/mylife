from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from numpy import array,savetxt,c_

f=open('/Users/dmelgar/Misc/IRIS_permanent.txt','r')
fout='/Users/dmelgar/Misc/IRIS_permanent_2015.txt'
tcompare=UTCDateTime('2015-12-31T23:59:59') #Stations that start after 2000 are EXCLDUDE
tcompare2=UTCDateTime('2015-01-01T00:00:00')
tcheck=UTCDateTime('1920-01-01T00:00:00.000000Z')
latout=[]
lonout=[]
k=0
while True:
    line=f.readline()
    if line=='':
        break
    if 'NETWORK' in line: #Ignore this line
        pass
    else:
        tstart=UTCDateTime(line.split('\t')[1])
        tend=UTCDateTime(line.split('\t')[2])
        if tstart>tcheck:
            if tstart<tcompare and tend>tcompare2: #Station starts before the endtime
                print str(tstart)+'   '+str(tend)
                k+=1
                latout.append(float(line.split('\t')[3]))
                lonout.append(float(line.split('\t')[4]))
print str(k)+' stations have data'
lonout=array(lonout)
latout=array(latout)
savetxt(fout,c_[lonout,latout],fmt='%12.4f\t%10.4f')
plt.figure()
plt.scatter(lonout,latout)
plt.title('Stations available by '+str(tcompare))
plt.show()