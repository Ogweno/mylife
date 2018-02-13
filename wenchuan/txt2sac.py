'''
D. Melgar 5/2015
Sample script to convert text data into SAC format usin ObsPy
'''

#Import functions I'm going to nead
from obspy import Stream,Trace,UTCDateTime
from numpy import genfromtxt 
from scipy.integrate import cumtrapz
from scipy.signal import filtfilt,butter
from matplotlib import pyplot as plt

#Paths and other global things
station_name='WCJYD'
path='/Users/dmelgar/Wenchuan2008/'
start_time=UTCDateTime('2008-05-12T06:28:25.068')
dt=0.005
fcorner=1./50 #High pass filter corner (1./50) will high pass filter above 50s period


#Read text data into variable
data_file=path+station_name+'.TXT'
acc=genfromtxt(data_file)
#Read into individual variabloes and convert to m/s/s
east_acc=acc[:,0]/100
north_acc=acc[:,1]/100
up_acc=acc[:,2]/100

#Create empty obspy stream objects (these will then be saved as SAC)
north=Stream(Trace())
east=Stream(Trace())
up=Stream(Trace())

#Put data in stream objects
north[0].data=north_acc
east[0].data=east_acc
up[0].data=up_acc

#Define the start time, sampling interval and station name
north[0].stats.starttime=start_time
north[0].stats.delta=dt
north[0].stats.station=station_name
east[0].stats.starttime=start_time
east[0].stats.delta=dt
east[0].stats.station=station_name
up[0].stats.starttime=start_time
up[0].stats.delta=dt
up[0].stats.station=station_name

#That's it! Now we apply post-processing, integrate to velocity and decimate to 5Hz
#
#Integrate
north[0].data=cumtrapz(north[0].data,north[0].times(),initial=0)
east[0].data=cumtrapz(east[0].data,east[0].times(),initial=0)
up[0].data=cumtrapz(up[0].data,up[0].times(),initial=0)

#High pass filter to remove baseline offsets
#This defines a high-pass filter function
def high_pass_filter(tr,fcorner,order):
    b, a = butter(order, fcorner,btype='highpass')
    y = filtfilt(b, a, tr)
    return y
def low_pass_filter(tr,fcorner,dt,order):
    from scipy.signal import filtfilt,butter
    fnyquist=1./(2*dt)
    print fnyquist
    Fc=fcorner/fnyquist
    print Fc
    b, a = butter(order, Fc,btype='lowpass')
    y = filtfilt(b, a, tr)
    return y
#Do the actual filtering
fnyquist=1./(2*east[0].stats.delta)
Fc=fcorner/fnyquist
filter_order=4
east[0].data=high_pass_filter(east[0].data,Fc,filter_order)
north[0].data=high_pass_filter(north[0].data,Fc,filter_order)
up[0].data=high_pass_filter(up[0].data,Fc,filter_order)

    
#This defines a decimating function with an anti-aliasing filter   
def stdecimate(st,factor,order):
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout

# We can't decimate from 200Hz to 5Hz in one step without introducing nuemrical error
# so we do it in three stages 200Hz -> 40Hz -> 10Hz -> 5Hz
north=stdecimate(north,5,13) ; north=stdecimate(north,4,13) ; north=stdecimate(north,2,13)
east=stdecimate(east,5,13) ; east=stdecimate(east,4,13) ; east=stdecimate(east,2,13)
up=stdecimate(up,5,13) ; u=stdecimate(up,4,13) ; u=stdecimate(up,2,13)

#Save to file as SAC format
north.write(path+station_name+'.n.vel',format='SAC')
east.write(path+station_name+'.e.vel',format='SAC')
up.write(path+station_name+'.u.vel',format='SAC')

#Done! Now make a plot to check that everything is fine
plt.figure()
plt.suptitle('Station '+north[0].stats.station)
plt.subplot(311)
plt.plot(north[0].times(),north[0].data)
plt.ylabel('North (m/s)')
plt.grid()

plt.subplot(312)
plt.plot(east[0].times(),east[0].data)
plt.ylabel('East (m/s)')
plt.grid()

plt.subplot(313)
plt.plot(up[0].times(),up[0].data)
plt.ylabel('up (m/s)')
plt.xlabel('Time (s)')
plt.grid()
plt.show()