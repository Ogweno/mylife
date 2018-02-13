def read_multi_segment_file(msf_file):
    
    from numpy import zeros
    
    #First count the number of segments and the number of elements in each segment
    Nsegments=0
    Num_elements=[]
    f=open(msf_file,'r')
    first_line=True
    while True:
        line=f.readline()
        if '>' in line:
            if first_line==False: #Append previous count
                Num_elements.append(Numel)
            first_line=False
            Nsegments+=1
            Numel=0
        else:
            Numel+=1
        if line=='': #End of file
            Num_elements.append(Numel-1)
            break
    f.close()
            
    #Now loop over segments and make an arra per segment adn append to list of arrays
    all_segments=[]    
    f=open(msf_file,'r')
    
    for ksegment in range(Nsegments):
        
        #First line in the segmetn is the stupid carrot
        line=f.readline()
        
        #Now read the next Num_elements[ksegment] lines
        lonlat=zeros((Num_elements[ksegment],2))
        for kelement in range(Num_elements[ksegment]):
            line=f.readline()
            lonlat[kelement,0]=float(line.split()[0])
            lonlat[kelement,1]=float(line.split()[1])
            
        #Done, append to list
        all_segments.append(lonlat)
    
    f.close()
    
    return all_segments
            
        
    