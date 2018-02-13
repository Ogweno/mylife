from mudpy import forward
from numpy import arange


#mudpy_file=u'/Users/dmelgar/Slip_inv/Chiapas/output/inverse_models/models/tsun.0012.inv'
#out_file=u'/Users/dmelgar/Slip_inv/Chiapas/output/inverse_models/models/tsun.0012.dtopo'
mudpy_file=u'/Users/dmelgar/Slip_inv/Chiapas/output/inverse_models/models/static_tsunami.0007.inv'
out_file=u'/Users/dmelgar/Slip_inv/Chiapas/output/inverse_models/models/static_tsunami.0007.dtopo'
x=arange(-97,-93,0.025)
y=arange(13,17,0.025)

forward.move_seafloor_okada(mudpy_file,out_file,x,y,refine_factor=None,mu=40e9,return_object=False)