from numpy import genfromtxt,zeros,where,std,nan

all_data_file='/Users/dmelgar/GPS_stuff/SoCal_xyz_24hr_dev1.log'
sta_solar_file='/Users/dmelgar/GPS_stuff/station_list/shallow_foundation_mast_solar_panel.sta'
sta_deepDB_file='/Users/dmelgar/GPS_stuff/station_list/deep_DB.sta'

sta_solar=genfromtxt(sta_solar_file,dtype='S')
sta_deepDB=genfromtxt(sta_deepDB_file,dtype='S')
all_stations=genfromtxt(all_data_file,usecols=0,dtype='S')
all_gps=genfromtxt(all_data_file,usecols=[2,3,4])

variance_solar=zeros(len(sta_solar))
variance_deepDB=zeros(len(sta_deepDB))
for k in range(len(sta_solar)):
    print k
    i=where(all_stations==sta_solar[k].upper())[0]
    if len(i)>0:
        print 'Has data'
        total=(all_gps[i,0]**2+all_gps[i,0]**2+all_gps[i,0]**2)**0.5
        variance_solar[k]=std(total/1e6)
    else:
        variance_solar[k]=nan
        
for k in range(len(sta_deepDB)):
    print k
    i=where(all_stations==sta_deepDB[k].upper())[0]
    if len(i)>0:
        print 'Has data'
        total=(all_gps[i,0]**2+all_gps[i,0]**2+all_gps[i,0]**2)**0.5
        variance_deepDB[k]=std(total/1e6)
    else:
        variance_deepDB[k]=nan
