from numpy import array,r_

fout=u'/Users/dmelgar/SSN/cataloges/201709.txt'
f=open(u'/Users/dmelgar/SSN/cataloges/201709.CAT')

k=0
while True:
    line=f.readline()
    if line=='':
        break
    if line[1:5]=='2017': #event
        if len(line.split())>11:
            print k
            lat=float(line[23:30])
            lon=float(line[30:38])
            z=float(line[38:43])
            time=line[1:20]
            if k==0:
                time_out=[time]
                lat_out=[lat]
                lon_out=[lon]
                z_out=[z]
            else:
                lat_out.append(lat)
                lon_out.append(lon)
                z_out.append(z)
                time_out.append(time)
            k+=1
f.close()

f=open(fout,'w')
for k in range(len(lon_out)):
    year=int(time_out[k][0:4])
    month=int(time_out[k][6])
    day=int(time_out[k][7:9])
    hour=int(time_out[k][10:12])
    minute=int(time_out[k][12:14])
    second=float(time_out[k][15:19])
    lon=float(lon_out[k])
    lat=float(lat_out[k])
    z=float(z_out[k])
    line='%4d %2d %2d %2d %2d %.1f\t%10.4f%10.4f\t%6.1f\n' %(year,month,day,hour,minute,second,lon,lat,z)
    f.write(line)
f.close()