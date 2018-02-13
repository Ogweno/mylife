from matplotlib import pyplot as plt
from glob import glob
from numpy import genfromtxt

path='/Users/dmelgar/Coquimbo2015/GPS/post/'
files=glob(path+'*.dat')

plt.figure()
for k in range(len(files)):
    day=genfromtxt(files[k],usecols=2)
    x=genfromtxt(files[k],usecols=3)
    y=genfromtxt(files[k],usecols=5)
    z=genfromtxt(files[k],usecols=7)
    plt.subplot(311)
    plt.plot(day,x-x[0])
    plt.ylabel(r'$\Delta x$')
    plt.subplot(312)
    plt.plot(day,y-y[0])
    plt.ylabel(r'$\Delta y$')
    plt.subplot(313)
    plt.plot(day,z-z[0])
    plt.ylabel(r'$\Delta z$')
    plt.xlabel('Day of year')
plt.show()
    