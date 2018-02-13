from mudpy import forward
from numpy import arange
#
#
#epi=[-83.520,17.469,13.0]
#fout='/Users/dmelgar/Honduras2018/fault/rapid.fault'
#forward.makefault(fout,254,65,50,3,3,epi,4,4,1.0)
#
#mudpy_file='/Users/dmelgar/Honduras2018/fault/rapid.rupt'
#out_file='/Users/dmelgar/Honduras2018/fault/rapid.dtopo'
#
#x=arange(-85,-81,0.025)
#y=arange(15.5,18.5,0.025)
#
#forward.move_seafloor_okada(mudpy_file,out_file,x,y,refine_factor=None,mu=30e9,return_object=False)


usgs_model='/Users/dmelgar/Honduras2018/fault/usgs.txt'
out_file_rupt='/Users/dmelgar/Honduras2018/fault/usgs.rupt'
Dx=5e3
Dy=3e3

x=arange(-87,-80,0.025)
y=arange(14.5,20.5,0.025)
forward.usgs2rupt(usgs_model,out_file_rupt,Dx,Dy)

out_file='/Users/dmelgar/Honduras2018/fault/usgs.dtopo'
forward.move_seafloor_okada(out_file_rupt,out_file,x,y,refine_factor=None,mu=30e9,return_object=False)

