from numpy import genfromtxt,where,c_,savetxt
out_path='/Users/dmelgar/code/GMT/coquimbo/contours/'

#Low-f parameters
scale=400
offset=0.08
name='0.02_0.5'
dotfile='/Users/dmelgar/Coquimbo2015/BP/20s/Chile2015.0.02_0.5Hz.Backprojection.timesteps.1s_step'

#high-f params
#scale=400
#offset=0.15
#name='0.5_2'
#dotfile='/Users/dmelgar/Coquimbo2015/BP/20s/Chile2015.0.5_2Hz.Backprojection.timesteps.1s_step'

dots=genfromtxt(dotfile)
dt=10

for k in range(10):
    i=where((dots[:,2]>=k*dt) & (dots[:,2]<(k+1)*dt))[0]
    out=c_[dots[i,4],dots[i,3],dots[i,2]-dots[i[0],2],(dots[i,5]/scale)+offset]
    outfile=out_path+'dots_'+name+'_'+str(dt*k)+'_'+str(dt*(k+1))+'s'
    savetxt(outfile,out,fmt='%.4f\t%.4f\t%.4f\t%.4f',header='lon,lat,power,time')