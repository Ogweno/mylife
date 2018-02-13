from numpy import genfromtxt,savetxt
from glob import glob

#contour_file=glob(u'/Users/dmelgar/Coquimbo2015/BP/10s/*0.5_2Hz*')
#contour_file=['/Users/dmelgar/Coquimbo2015/BP/10s/Chile2015.0.02_0.5Hz.Backprojection.contour.total']
#contour_file=['/Users/dmelgar/Coquimbo2015/BP/10s/Chile2015.0.5_2Hz.Backprojection.contour.total']
#contour_file=glob(u'/Users/dmelgar/Coquimbo2015/BP/20s/*0.5_2Hz*')
contour_file=glob(u'/Users/dmelgar/Coquimbo2015/BP/20s/*0.02_0.5Hz*')
out_path='/Users/dmelgar/code/GMT/coquimbo/contours/'
name='0.02_0.5'

for kc in range(len(contour_file)):
    data=genfromtxt(contour_file[kc])
    k=0
    while True:
        contour=data[k,0]
        num_points=data[k,1]
        print contour
        k+=1
        file_out=out_path+name+'_'+contour_file[kc].split('.')[-1]+'_'+str(int(contour))
        if k==1:
            out=data[k:k+num_points,:]
        else:
            out=data[k:k+num_points,:]
        savetxt(file_out,out,fmt='%.6f\t%.6f')
        k=k+num_points
        if k>=len(data):
            break
    