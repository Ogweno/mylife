from glob import glob
from numpy import float64

path='/Users/dmelgar/Coquimbo2015/EEW/2015/'
dataless_path='/Users/dmelgar/Coquimbo2015/EEW/dataless/'
outfile='/Users/dmelgar/Coquimbo2015/EEW/stations.txt'
fout=open(outfile,'w')
fout.write('# NET,STA,CHAN,LON,LAT,GAIN (counts->m/s or counts->m/s/s)\n')

networks=glob(path+'*')
for net in networks:
    stations=glob(net+'/*')
    for sta in stations:
        channels=glob(sta+'/*')
        for chan in channels:
            #Get sncl
            sncl_resp=net.split('/')[-1]+'.'+sta.split('/')[-1]+'..'+chan.split('/')[-1].split('.')[0]
            sncl_pz=net.split('/')[-1]+'_'+sta.split('/')[-1]+'_'+chan.split('/')[-1].split('.')[0]
            print sncl_resp
            #Now go look for response file
            resp_file=glob(dataless_path+'*'+sncl_resp+'*')
            pz_file=glob(dataless_path+'*'+sncl_pz+'*')
            if len(resp_file)>0:  #Parse it
                # Get latitude and longitude
                f=open(pz_file[0],'r')
                while True:
                    line=f.readline()
                    if 'latitude' in line.lower():
                        lat=float(line.split()[-1])
                    if 'longitude' in line.lower():
                        lon=float(line.split()[-1])
                    if 'elevation' in line.lower():
                        elev=float(line.split()[-1])
                    if not line:
                        break
                f.close()
                # Get simple gain
                f=open(resp_file[0],'r')
                found_sensitivity=False
                stage_0=False
                infinite=0
                while True:
                    line=f.readline()
                    if 'Channel Sensitivity' in line:
                        found_sensitivity=True
                    if found_sensitivity:
                        if 'Stage sequence number:                 0' in line:
                            stage_0=True
                        if stage_0==True:
                            if 'Sensitivity' in line:
                                gain=float64(line.split()[-1])
                                break
                    if infinite>5000:
                        print 'Infinite loop'
                        break
                    infinite+=1
                f.close()
                #Make output to table file
                out = '%2s\t%4s\t%3s\t%10.6f\t%10.6f\t%6i\t%12i\n' %(sncl_pz.split('_')[0],sncl_pz.split('_')[1],sncl_pz.split('_')[2],lon,lat,int(elev),gain)
                fout.write(out)
            else:
                print 'No RESP file found for '+sncl_resp
            # Write to gain table file
fout.close()
            
                
        