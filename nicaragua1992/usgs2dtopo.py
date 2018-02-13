from mudpy import forward
from numpy import arange

Dx=12e3
Dy=8.10e3
mu=10e9

usgs_model='/Users/dmelgar/Nicaragua1992/p0005ddn.param'
mudpy_file='/Users/dmelgar/Nicaragua1992/nicaragua.rupt'
dtopo_file='/Users/dmelgar/Nicaragua1992/nicaragua.dtopo'

forward.usgs2rupt(usgs_model,mudpy_file,Dx,Dy)


x=arange(-91,-84,0.02)
y=arange(8,14,0.02)

forward.move_seafloor_okada(mudpy_file,dtopo_file,x,y,mu=mu)