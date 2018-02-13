from numpy import zeros,savetxt
import re

kml='/Users/dmelgar/KMLs/japan/JPN_adm1.kml'
f=open(kml,'r')
prefecture=''
while True:
    if prefecture=='Akita':
        a=0
    line=f.readline()
    if 'name="NAME_1"' in line:
        prefecture=line.split('>')[-2].split('<')[0]
        print prefecture
    if '<MultiGeometry>' in line:
        #Parse line
        parser=re.finditer('<coordinates>(.+?)</coordinates>', line)
        coords_out=[]
        for coords in parser:
            coords_out.append(coords.group(1))
        #Open file
        fout='/Users/dmelgar/KMLs/japan/prefectures/JPN_'+prefecture+'.txt'
        fpref=open(fout,'w')
        #Now split each one and put in multi segment file GMT can understand
        for k in range(len(coords_out)):
            fpref.write('>\n')
            coord_string=coords_out[k].split(' ')
            for kwrite in range(len(coord_string)):
                line_out=coord_string[kwrite].split(',')[0]+' '+coord_string[kwrite].split(',')[1]+'\n'
                fpref.write(line_out)
        fpref.close()
    if line=='':
        break
            
            

            

