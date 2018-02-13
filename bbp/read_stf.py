
from numpy import array,r_,zeros

#f=open('/Users/dmelgar/code/BBP/bbp/bbp_data/indata/7301769/northridge_eq_gp.srf','r')
f=open(u'/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/65/test_eew_m6.5.srf')
 
begin_read=False
line=f.readline()
kpoint=0
while True:
    if begin_read==False:
        line=f.readline()
    if 'POINTS' in line:
        points=int(line.split()[1])
        begin_read=True
        xyz=zeros((points,3))
        slip=zeros(points)
        tinit=zeros(points)
    if begin_read:
        #Read infor about this point source location, type, etc
        line1=f.readline()
        if line1 == '':
            break
        #Assign latlon points
        xyz[kpoint,0]=float(line1.split()[0])
        xyz[kpoint,1]=float(line1.split()[1])
        xyz[kpoint,2]=float(line1.split()[2])
        tinit[kpoint]=float(line1.split()[6])
        line2=f.readline()
        slip[kpoint]=float(line2.split()[1])
        kpoint+=1
        #Number of STF points
        Nstf=int(line2.split()[2])
        #Format only allows 6 pts epr line, so how many lines do we need to read?
        if Nstf % 6==0:
            Nread=Nstf/6
        else:
            Nread=Nstf/6+1
        #Read stf for this point source
        stf=array([])
        for kread in range(Nread):
            line=f.readline()
            stf=r_[stf,array(line.split()).astype('float')]

     
f.close()

