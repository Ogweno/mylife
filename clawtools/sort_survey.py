from numpy import argsort,genfromtxt,savetxt

s=genfromtxt('/Users/dmelgar/Maule2010/tsunami/survey.txt')
i=argsort(s[:,1])
lat=s[i,1]
lon=s[i,2]
height=s[i,3]
s[:,1]=lat
s[:,2]=360+lon
s[:,3]=height
savetxt('/Users/dmelgar/Maule2010/tsunami/survey_sort.txt',s,fmt='%d\t%.6f\t%.6f\t%.6f')