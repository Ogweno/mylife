from mudpy.forward import move_seafloor_okada
from numpy import arange

mudpy_file='/Users/dmelgar/Slip_inv/Wharton_tsun/forward_models/model6c.rupt'
out_file='/Users/dmelgar/Slip_inv/Wharton_tsun/forward_models/model6c.dtopo'
x=arange(89,95,0.01)
y=arange(0.5,4.0,0.01)
refine_factor=None
mu=40e9


fault=move_seafloor_okada(mudpy_file,out_file,x,y,refine_factor=refine_factor,mu=mu)
