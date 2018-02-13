import requests
from numpy import genfromtxt


stations_noaa=genfromtxt('/Users/dmelgar/tidegauge_noise/station_info/stations.txt',usecols=0,dtype='S')
stations_txt=genfromtxt('/Users/dmelgar/tidegauge_noise/station_info/stations.txt',usecols=1,dtype='S')

path_out='/Users/dmelgar/tidegauge_noise/data/'

years=['2017']
months=['01','02','03','04','05','06','07','08','09','10','11','12']
days=['31','28','31','30','31','30','31','31','30','31','30','31']


#variable='tide'
#variable='pressure'
variable='temperature'

for ksta in range(len(stations_noaa)):
    
    print 'Working on station '+stations_txt[ksta]
    
    for kmonth in range(len(months)):
        
        print '... month '+months[kmonth]
        
        #Form the web request url
        t1=years[0]+months[kmonth]+'01'
        t2=years[0]+months[kmonth]+days[kmonth]
        sta=stations_noaa[ksta]
        
        if variable=='tide':
            url='https://tidesandcurrents.noaa.gov/api/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&begin_date='+t1+'&end_date='+t2+'&datum=MLLW&station='+sta+'&time_zone=GMT&units=metric&format=csv'
            suffix=''
        elif variable=='temperature':
            url='https://tidesandcurrents.noaa.gov/api/datagetter?product=air_temperature&application=NOS.COOPS.TAC.WL&begin_date='+t1+'&end_date='+t2+'&station='+sta+'&time_zone=GMT&units=metric&format=csv'
            suffix='temp'
        elif variable=='pressure':
            url='https://tidesandcurrents.noaa.gov/api/datagetter?product=air_pressure&application=NOS.COOPS.TAC.WL&begin_date='+t1+'&end_date='+t2+'&station='+sta+'&time_zone=GMT&units=metric&format=csv'
            suffix='pres'
        else:
            print 'ERROR: don;t know what that variable is.';
            break
            
            
        r=requests.get(url)
        data=r.text
        
        #Replace ',' with \t
        data=data.replace(',','\t')
        
        #Add header line hash
        data=data.replace('Date','#Date')
        
        #Output folder and write
        fout=open(path_out+stations_txt[ksta]+'/csv/'+stations_txt[ksta]+'_'+suffix+'_'+months[kmonth]+'_'+years[0]+'.csv','w')
        fout.write(data)
        fout.close()
        
        
        
        
    