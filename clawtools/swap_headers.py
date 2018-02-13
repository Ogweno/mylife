#invert .topo headers

from glob import glob
from fileinput import FileInput
from os import rename

#files=glob(u'/Users/dmelgar/DEMs/iquique/survey/*.topo')
#tmpname='/Users/dmelgar/DEMs/iquique/survey/tmp'
files=glob(u'/Users/dmelgar/DEMs/iquique/gauges/*.topo')
tmpname='/Users/dmelgar/DEMs/iquique/gauges/tmp'

for k in range(len(files)):
    linecount=0
    f1=open(files[k],'r')
    f2=open(tmpname,'w')
    while True:
        line=f1.readline()
        if linecount==0:
            ncols=line.split()[1]
            lineout=ncols+"\tmx\n"
        elif linecount==1:
            nrows=line.split()[1]
            lineout=nrows+"\tmy\n"
        elif linecount==2:
            xl=line.split()[1]
            if xl>0:
                xl=repr(float(xl)-360)
            lineout=xl+"\txllcenter\n"
        elif linecount==3:
            yl=line.split()[1]
            lineout=yl+"\tyllcenter\n"
        elif linecount==4:
            cell=line.split()[1]
            lineout=cell+"\tcellsize\n"
        elif linecount==5:
            nan=line.split()[1]
            lineout=nan+"\tnodata_value\n"
        elif linecount>5:
            lineout=line
        if line=='':
            break
        f2.write(lineout)
        linecount+=1
    f1.close()
    f2.close()
    rename(tmpname,files[k])