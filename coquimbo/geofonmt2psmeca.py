from glob import glob

files=glob('/Users/dmelgar/Coquimbo2015/afters/mt-geofon/2015*.txt')
fout=u'/Users/dmelgar/code/GMT/coquimbo/MTs.txt'
out=''
for k in range(len(files)):
    f=open(files[k],'r')
    print files[k]
    while True:
        l=f.readline()
        if 'Epicenter' in l:
            lat=float(l.split(' ')[1])
            lon=float(l.split(' ')[2])
        if 'Depth' in l:
            z=float(l.split()[1])
        if 'Scale' in l:
            scale=int(l.split(' ')[5].split('**')[1])
        if 'Mrr' in l:
            mrr=float(l.split("=")[1].split()[0])
            mtt=float(l.split("=")[2])
        if 'Mpp' in l:
            mff=float(l.split("=")[1].split()[0])
            mrt=float(l.split("=")[2])
        if 'Mrp' in l:
            mrf=float(l.split("=")[1].split()[0])
            mtf=float(l.split("=")[2])
            break
    psmeca='%.2f\t%.2f\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\n' %(lon,lat,z,mrr,mtt,mff,mrt,mrf,mtf,scale)
    out=out+psmeca
    f.close()
f2=open(fout,'w')
f2.write(out)
f2.close()