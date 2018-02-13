'''
Make a fault model. First we will make a bialteral fault model. That means it has 
the same number of subfaults on both sides of the hypocenter. But the Wenchuan
earthquake only ruptured to the north of the hypocenter so then we will cut the
model to only include more faults north of the hypocenter than to the south
'''


from mudpy.forward import makefault
from numpy import array,genfromtxt,arange,r_,savetxt
from matplotlib import pyplot as plt 

#Basic geoemtry
strike=228
dip=38
#where to save the .fault file
fout='/Users/dmelgar/Slip_inv/wenc_2008/data/model_info/wenc_bilateral.fault'
fout2='/Users/dmelgar/Slip_inv/wenc_2008/data/model_info/wenc_unilateral.fault'
#Number of subfaults along strike
nstrike=74
#Number of subfaults ABOVE the hypocenter
num_updip=2
#Number of subfaults BELOW the hypocenter
num_downdip=5
#Define size of subfaults
dx_strike=7 #In km
dx_dip=7 #In km
#
epicenter=array([103.36,31.0,12])  #I used the USGS one
#Rise time
rise_time=3.0
#Make the fault
makefault(strike,dip,nstrike,dx_dip,dx_strike,epicenter,num_updip,num_downdip,rise_time,fout)

#Read the model you just made so you can make a plot and some editting
fault=genfromtxt(fout)
#Make a quick plot so you can see the model
plt.figure()
plt.scatter(fault[:,1],fault[:,2])
plt.scatter(epicenter[0],epicenter[1],c='r',s=80)
plt.title('Bilateral model')
plt.show()

# Ok now make unialteral model, we will keep all subfaults north of the hypocenter
# but we will only keep 2 columns of faults to the south. This is determined
# by looking at the plot of the subfaults

ndip=num_updip+num_downdip+1 #Total number of rows of subfaults
good_faults=arange(0,nstrike-18) #We want tot hrow away the first 18 columns of subfaults
keep_faults=good_faults.copy() #Initialize
for k in range(ndip-1):
    keep_faults=r_[keep_faults,good_faults+nstrike*(k+1)]

#Ok now throw away faults you don't want
fault=fault[keep_faults,:]
#Renumber the faults so the numbers are consecutive
fault[:,0]=arange(1,len(fault)+1)
#Save to file
savetxt(fout2,fault,fmt='%i\t%.6f\t%.6f\t%.3f\t%i\t%i\t%.1f\t%.1f\t%.2f\t%.2f')

#Make a quick plot so you can see the model
plt.figure()
plt.scatter(fault[:,1],fault[:,2])
plt.scatter(epicenter[0],epicenter[1],c='r',s=80)
plt.title('Unilateral model')
plt.show()
