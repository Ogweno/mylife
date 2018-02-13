from mudpy import insartools
from mudpy import forward



#process the quadtree file

#home='/Users/dmelgar/Slip_inv/'
#project_name='Ixtepec'
#file_in=u'/Users/dmelgar/Ixtepec2017/InSAR/T107_Ascending/ascending_quadtree.los'
#gf_list=home+project_name+'/data/station_info/ascending.gflist'
#prefix='AS'
#
#insartools.quadtree2mudpy(home,project_name,file_in,gf_list,prefix)


# making the fault geometry

fault_file='/Users/dmelgar/Slip_inv/Ixtepec/data/model_info/usgs_NP1.fault'
strike=79
dip=55
num_columns=20
dx_strike=2.0
dx_dip=2.0
hypocenter=[-95.168597, 16.494261,10.0]
number_updip=5
number_downdip=5
rise_time=1.0

forward.makefault(fault_file,strike,dip,num_columns,dx_strike,
                    dx_dip,hypocenter,number_updip,number_downdip,rise_time)

