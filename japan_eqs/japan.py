from numpy import genfromtxt,zeros
from datetime import datetime
from matplotlib import pyplot as plt

t=genfromtxt('/Users/dmelgar/code/python/japan_eqs/catalogue.txt',usecols=0,dtype='S')
year=zeros(len(t))
for k in range(len(t)):
    temp=datetime.strptime(t[k],'%Y/%m/%d,%H:%M:%S.%f')
    year[k]=temp.year+(temp.month-1)*(1./12)+temp.day*(1./365)
M=genfromtxt('/Users/dmelgar/code/python/japan_eqs/catalogue.txt',usecols=11)


plt.close("all")
plt.figure()
plt.subplot(121)
plt.stem(year,M)
plt.ylim([6,9.5])
plt.xlim([2000,2015])
plt.xlabel(r'Decimal Year')
plt.ylabel(r'Moment Mangitude')
plt.suptitle(r'NIED Catalogue, 82 events with $M_w$>6.5',fontsize=18)
plt.grid()
plt.subplot(122)
n, bins, patches = plt.hist(M, 5, facecolor='green')
plt.xlabel('Moment Magnitude')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()
