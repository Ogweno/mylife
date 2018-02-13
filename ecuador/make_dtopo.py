from mudpy.forward import move_seafloor_okada
from numpy import arange

#mudpy_file='/Users/dmelgar/Slip_inv/Ecuador_insar/output/inverse_models/models/sm_insar_5win_vr3.4.0014.noshallow.inv.total'
#out_file='/Users/dmelgar/Ecuador2016/tsunami/dtopo/sm_insar_5win_vr3.4.0014.noshallow.dtopo'
mudpy_file='/Users/dmelgar/Slip_inv/Ecuador_insar/output/inverse_models/models/sm_insar_5win_vr3.4.0014.inv.total'
out_file='/Users/dmelgar/Ecuador2016/tsunami/dtopo/sm_insar_5win_vr3.4.0014.dtopo'
x=arange(-82,-79,0.005)
y=arange(-1.5,1.5,0.005)
refine_factor=None
mu=40e9


fault=move_seafloor_okada(mudpy_file,out_file,x,y,refine_factor=refine_factor,mu=mu,return_object=True)
