from glob import glob

path_all=u'/Users/dmelgar/Coquimbo2015/strong_motion/raw/'
stations=glob(path_all+'*HNN*')

for k in range(len(stations)):
    sta=stations[k].split('.')[7]
    network=stations[k].split('.')[6]
    print sta+' , '+network
    if network=='C':
        resp=u'/Users/dmelgar/Coquimbo2015/strong_motion/red_C/red_C.dataless_'+sta
    else:
        resp=u'/Users/dmelgar/Coquimbo2015/strong_motion/red_C1/red_C1.dataless_'+sta
        