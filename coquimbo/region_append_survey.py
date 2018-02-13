from numpy import genfromtxt
#coast=genfromtxt(u'/Users/dmelgar/Maule2010/tsunami/maule_3_coast_fg.txt')
#f=open(u'/Users/dmelgar/Maule2010/tsunami/gauge_region.txt','w')
#dl=6*3./3600
#for k in range(len(coast)):
#    #rundata.regiondata.regions.append([5, 5, 0, 1e10, -70.967500000000, -70.817500000000, -18.071666666667, -18.000000000000])
#    out= 'rundata.regiondata.regions.append([5, 5, 0, 1e10,'+str(coast[k,0]-dl)+','+str(coast[k,0]+dl)+','+str(coast[k,1]-dl)+','+str(coast[k,1]+dl)+'])\n'
#    f.write(out)
#f.close()

#survey=genfromtxt(u'/Users/dmelgar/Tsunamis/coquimbo_3_gps_sm_tg/fgrid_coast.txt',skip_header=7)
survey=genfromtxt(u'/Users/dmelgar/Tsunamis/coquimbo_30_gps_sm_tg/fgrid_coast.txt',skip_header=7)
dl=5*30./3600
f=open(u'/Users/dmelgar/Coquimbo2015/tsunami/coast_append_region.txt','w')
for k in range(len(survey)):
    #rundata.regiondata.regions.append([5, 5, 0, 1e10, -70.967500000000, -70.817500000000, -18.071666666667, -18.000000000000])
    out= 'rundata.regiondata.regions.append([3, 3, 0, 1e10,'+str(survey[k,0]-dl)+','+str(survey[k,0]+dl)+','+str(survey[k,1]-dl)+','+str(survey[k,1]+dl)+'])\n'
    f.write(out)
f.close()
