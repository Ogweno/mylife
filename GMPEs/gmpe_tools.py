'''
Assorted ground motion prediction equation tools
Diego Melgar, May, 2016
'''


def PGAr_calc(M, Rjb, U, SS, RS, NS,italy=False):
    '''
    Calculate reference PGA
    '''
    
    from numpy import log,exp,ones,where
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    
    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    
    Mharray=Mh*ones(Rjb.shape)
    Uarray=U.copy()
    SSarray=SS.copy()
    NSarray=NS.copy()
    RSarray=RS.copy()
    harray=h*ones(Rjb.shape)
    c1array=c1*ones(Rjb.shape)
    c2array=c2*ones(Rjb.shape)
    Mrefarray=Mref*ones(Rjb.shape)
    c3array=c3*ones(Rjb.shape)
    Dc3array=Dc3*ones(Rjb.shape)
    Dc3jpitarray=Dc3jpit*ones(Rjb.shape)
    Rrefarray=Rref*ones(Rjb.shape)
    
    fm=ones(Rjb.shape)
    i=where((M<Mharray)==True)[0]
    #if M <= Mh:
    fm[i] = e0*Uarray[i] + e1*SSarray[i] + e2*NSarray[i] + e3*RSarray[i] + e4*(M[i] - Mharray[i]) + e5*(M[i] - Mharray[i])**2
    #else:
    i=where((M<Mharray)==False)[0]
    fm[i] = e0*Uarray[i] + e1*SSarray[i] + e2*NSarray[i] + e3*RSarray[i] + e6*(M[i] - Mharray[i])
        
    
    R = (Rjb**2 + harray**2)**0.5
    
    #region term
    #CA
    fp = (c1array + c2array * (M - Mrefarray)) * log(R / Rrefarray) + (c3array + Dc3array) * (R - Rrefarray)
    #ITLAY
    if italy==True:
        fp = (c1array + c2array * (M - Mrefarray)) * log(R / Rrefarray) + (c3array + Dc3jpitarray) * (R - Rrefarray)
    
    #Calculate PGAr
    PGAr = exp(fm + fp)
    
    return PGAr




def PGAr_calc_one_station(M, Rjb, U, RS, NS):
    '''
    Calculate reference PGA
    '''
    
    from numpy import log,exp,ones
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    
    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    
    if NS == 0 and RS == 0 and U == 0:
        SS = 1
    else:
        SS = 0
    

    if M <= Mh:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e4*(M - Mh) + e5*(M - Mh)**2
    else:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e6*(M - Mh)
    
    R = (Rjb**2 + h**2)**0.5
    
    #region term
    fp = (c1 + c2 * (M - Mref)) * log(R / Rref) + (c3 + Dc3) * (R - Rref)
    
    #Calculate PGAr
    PGAr = exp(fm + fp)
    
    return PGAr


        

def bssa14(M, Rjb, Vs30,U=None,SS=None,RS=None,NS=None,Z1=None,intensity_measure='PGA',italy=False):
    '''
    Calculate ground motion intensity using the BSSA14 GMPE
    
    Parameters:
        M - Moment magnitude
        Rjb - Distance to surface projection of fault in km
        U - is 1 if unspecified faulting style
        RS - is 1 if reverse faulting
        NS - is 1 if normal fault
        Vs30 - Vs30 in m/s
        Z1 - Depth to Vs=1km/s, if unknown use Z1=None
        
    Returns:
        Y - the desired ground motion intensity, PGA in g or PGV in cm/s
        
    Notes: For strike slip faulting (default) set U=NS=RS=0
    '''
    
    from numpy import log,exp,sqrt,array,ones,where,zeros
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    
    #Convert input to floats
    Vs30=Vs30.astype(float)
    Rjb=Rjb.astype(float)
    M=M.astype(float)

    
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    elif intensity_measure.upper()=='SA0.3':
        coefficients=[1.2217,1.2401,1.0246,1.2653,0.95676,-0.1959,-0.092855,6.14,-1.09480,0.13388,-0.00548,4.5,
                        1,4.93,0.000000,0.002200,-0.003300,-0.8417,1308.47,760,0,0.1,-0.2191,-0.00670,-9.9,-9.9,
                        103.150,268.590,0.138,0.050,225,300,0.6750,0.5610,0.3630,0.2290]
    elif intensity_measure.upper()=='SA1.0':
        coefficients=[0.3932,0.4218,0.207,0.4124,1.5004,-0.18983,0.17895,6.2,-1.19300,0.10248,-0.00121,4.5,1,5.74,
                        0.000000,0.002920,-0.002090,-1.0500,1109.95,760,0,0.1,-0.1052,-0.00844,0.367,0.208,116.390,
                        270.000,0.098,0.020,225,300,0.5530,0.6250,0.4980,0.2980]
    elif intensity_measure.upper()=='SA3.0':
        coefficients=[-1.1898,-1.142,-1.23,-1.2664,2.1323,-0.04332,0.62694,6.2,-1.21790,0.09764,0.00000,4.5,1,6.93,
                        0.000000,0.002620,-0.001190,-1.0112,922.43,760,0,0.1,-0.0136,-0.00183,1.135,0.516,130.360,
                        195.000,0.088,0.000,225,300,0.5340,0.6190,0.5370,0.3440]
    else:
        print 'ERROR: Unknown intensity measure'
        #return

    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    C = coefficients[17]
    Vc = coefficients[18]
    Vref = coefficients[19]
    f1 = coefficients[20]
    f3 = coefficients[21]
    f4 = coefficients[22]
    f5 = coefficients[23]
    f6 = coefficients[24]
    f7 = coefficients[25]
    R1 = coefficients[26]
    R2 = coefficients[27]
    Dfr = coefficients[28]
    Dfv = coefficients[29]
    V1 = coefficients[30]
    V2 = coefficients[31]
    phi1 = coefficients[32]
    phi2 = coefficients[33]
    tau1 = coefficients[34]
    tau2 = coefficients[35]


    Uarray=U.copy()
    NSarray=NS.copy()
    SSarray=SS.copy()
    RSarray=RS.copy()

    # Hinge magnitude term and conversion to arrays
    Mharray=Mh*ones(Rjb.shape)

    fm=zeros(Rjb.shape)
    
    i=where((M<=Mharray)==True)[0]
    #if M <= Mh:
    fm[i] = e0*Uarray[i] + e1*SSarray[i] + e2*NSarray[i] + e3*RSarray[i] + e4*(M[i] - Mharray[i]) + e5*(M[i] - Mharray[i])**2
    #else:
    i=where((M<=Mharray)==False)[0]
    fm[i] = e0*Uarray[i] + e1*SSarray[i] + e2*NSarray[i] + e3*RSarray[i] + e6*(M[i] - Mharray[i])
    
    #Disance term
    harray=h*ones(Rjb.shape)
    R = (Rjb**2 + harray**2)**0.5

    # Region term
    c1array=c1*ones(Rjb.shape)
    c2array=c2*ones(Rjb.shape)
    Mrefarray=Mref*ones(Rjb.shape)
    Rrefarray=Rref*ones(Rjb.shape)
    c3array=c3*ones(Rjb.shape)
    Dc3array=Dc3*ones(Rjb.shape)
    Dc3jpitarray=Dc3jpit*ones(Rjb.shape)
    if italy==True:
        fp = (c1array + c2array*(M - Mrefarray))*log(R/Rrefarray) + (c3array + Dc3jpitarray)*(R - Rrefarray)
    else:
        fp = (c1array + c2array*(M - Mrefarray))*log(R/Rrefarray) + (c3array + Dc3array)*(R - Rrefarray)

    #Linear Site Term
    Carray=C*ones(Rjb.shape)
    Vcarray=Vc*ones(Rjb.shape)
    Vrefarray=Vref*ones(Rjb.shape)
    i_vs30=(Vs30 <= Vc)
    i=where(i_vs30==True)[0]
    if len(i)>0:
        flin = Carray*log(Vs30[i]/Vrefarray[i])
    i=where(i_vs30==False)[0]
    if len(i)>0:
        flin = Carray*log(Vcarray[i]/Vrefarray[i])

    #Nonlinear Site Term    
    minVarray=Vs30.copy()
    i=where(Vs30>760)[0]
    if len(i)>0:
        minVarray[i_vs30] = 760

    #Combine terms
    PGAr=PGAr_calc(M, Rjb, U,SS, RS, NS)
    
    f5array=f5*ones(Rjb.shape)
    f1array=f1*ones(Rjb.shape)
    f3array=f3*ones(Rjb.shape)
    f4array=f4*ones(Rjb.shape)
    
    f2array = f4array*((exp(f5array*(minVarray - 360))) - exp(f5array*(760 - 360)))
    fnl = f1array + f2array*log((PGAr + f3array)/f3array)
    fnl = f1array + (f4array*((exp(f5array*(minVarray - 360)))-exp(f5*(760 - 360))))*log((PGAr + f3array)/f3array)

    #Basin Depth Term
    mz1 = exp(-7.15/4*log((Vs30**4 + 570.94**4)/(1360**4 + 570.94**4)))/1000

    #Final correction
    if Z1 == None:
        dz1 = zeros(Vs30.shape)
    else:
        dz1 = Z1 - mz1

    fz1 = zeros(Vs30.shape)
    
    #elif dz1 <= f7/f6:
    #    fz1 = f6*dz1
    #elif dz1 > f7/f6:
    #    fz1 = f7
    #else:
    #    fz1 = 0

    if Z1 == None:
        fz1 = zeros(Vs30.shape)
    else:
        fz1 = fz1

    #Site Term
    fs = flin + fnl #in ln units

    #Model Prediction in ln units
    
    Y = exp(fm + fp + fs + fz1)
    
    #Stdev
    sigma=bssa14_stdev(M,Rjb,Vs30,intensity_measure=intensity_measure)
    
    return Y,sigma

def bssa14_one_station(M, Rjb, Vs30, U=0, RS=0, NS=0, SS=0,Z1=None,intensity_measure='PGA'):
    '''
    Calculate ground motion intensity using the BSSA14 GMPE
    
    Parameters:
        M - Moment magnitude
        Rjb - Distance to surface projection of fault in km
        U - is 1 if unspecified faulting style
        RS - is 1 if reverse faulting
        NS - is 1 if normal fault
        Vs30 - Vs30 in m/s
        Z1 - Depth to Vs=1km/s, if unknown use Z1=None
        
    Returns:
        Y - the desired ground motion intensity, PGA in g or PGV in cm/s
        
    Notes: For strike slip faulting (default) set U=NS=RS=0
    '''
    
    from numpy import log,exp,sqrt
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    
    #Convert input to floats
    Vs30=float(Vs30)
    Rjb=float(Rjb)
    M=float(M)
    
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    else:
        print 'ERROR: Unknown intensity measure'
        #return

    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    C = coefficients[17]
    Vc = coefficients[18]
    Vref = coefficients[19]
    f1 = coefficients[20]
    f3 = coefficients[21]
    f4 = coefficients[22]
    f5 = coefficients[23]
    f6 = coefficients[24]
    f7 = coefficients[25]
    R1 = coefficients[26]
    R2 = coefficients[27]
    Dfr = coefficients[28]
    Dfv = coefficients[29]
    V1 = coefficients[30]
    V2 = coefficients[31]
    phi1 = coefficients[32]
    phi2 = coefficients[33]
    tau1 = coefficients[34]
    tau2 = coefficients[35]

    # Magnitude Scaling Term
    if NS == 0 and RS == 0 and U == 0:
        SS = 1
    else:
        SS = 0

    # Hinge magnitude term
    if M <= Mh:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e4*(M - Mh) + e5*(M - Mh)**2
    else:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e6*(M - Mh)  
    
    #Disance term
    R = (Rjb**2 + h**2)**0.5

    # Region term
    fp = (c1 + c2*(M - Mref))*log(R/Rref) + (c3 + Dc3)*(R - Rref)

    #Linear Site Term
    if Vs30 <= Vc:
        flin = C*log(Vs30/Vref)
    else:
        flin = C*log(Vc / Vref)

    #Nonlinear Site Term
    if Vs30 < 760:
        minV = Vs30
    else:
        minV = 760

    #Combine terms
    PGAr=PGAr_calc_one_station(M, Rjb, U, RS, NS)
    
    f2 = f4*((exp(f5*(minV - 360))) - exp(f5*(760 - 360)))
    fnl = f1 + f2*log((PGAr + f3)/f3)
    fnl = f1 + (f4*((exp(f5*(minV - 360)))-exp(f5*(760 - 360))))*log((PGAr + f3)/f3)

    #Basin Depth Term
    mz1 = exp(-7.15/4*log((Vs30**4 + 570.94**4)/(1360**4 + 570.94**4)))/1000

    #Final correction
    if Z1 == None:
        dz1 = 0
    else:
        dz1 = Z1 - mz1

    fz1 = 0
    
    #elif dz1 <= f7/f6:
    #    fz1 = f6*dz1
    #elif dz1 > f7/f6:
    #    fz1 = f7
    #else:
    #    fz1 = 0

    if Z1 == None:
        fz1 = 0
    else:
        fz1 = fz1

    #Site Term
    fs = flin + fnl #in ln units

    #Model Prediction in ln units
    
    Y = exp(fm + fp + fs + fz1)
    
    #Standard deviation
    sigma=bssa14_stdev_one_station(M,Rjb,Vs30,intensity_measure=intensity_measure)
    
    return Y,sigma


def bssa14_stdev_one_station(M,Rjb,Vs30,intensity_measure='PGA'):
    '''
    Get GMPE standar deviation
    '''
    from numpy import log,exp,sqrt,array,ones,where,zeros
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    
    #Convert input to floats
    Vs30=float(Vs30)
    Rjb=float(Rjb)
    M=float(M)
    
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    else:
        print 'ERROR: Unknown intensity measure'
        #return

    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    C = coefficients[17]
    Vc = coefficients[18]
    Vref = coefficients[19]
    f1 = coefficients[20]
    f3 = coefficients[21]
    f4 = coefficients[22]
    f5 = coefficients[23]
    f6 = coefficients[24]
    f7 = coefficients[25]
    R1 = coefficients[26]
    R2 = coefficients[27]
    Dfr = coefficients[28]
    Dfv = coefficients[29]
    V1 = coefficients[30]
    V2 = coefficients[31]
    phi1 = coefficients[32]
    phi2 = coefficients[33]
    tau1 = coefficients[34]
    tau2 = coefficients[35]


    if M <= 4.5:
        tauM = tau1
    elif M > 4.5 and M < 5.5:
        tauM = tau1 + (tau2 - tau1) * (M - 4.5)
    else:
        tauM = tau2
    
    if M <= 4.5:
        phiM = phi1
    elif M > 4.5 and M < 5.5:
        phiM = phi1 + (phi2 - phi1) * (M - 4.5)
    else:
        phiM = phi2

    
    if Rjb <= R1:
        phiMR = phiM
    elif Rjb > R1 and Rjb <= R2:
        phiMR = phiM + Dfr * (log(Rjb / R1) / (log(R2 / R1)))
    else:
        phiMR = phiM + Dfr
    
    if Vs30 >= V2:
        phi = phiMR
    elif Vs30 >= V1 and Vs30 <= V2:
        phi = phiMR - Dfv * (log(V2 / Vs30) / (log(V2 / V1)))
    else:
        phi = phiMR - Dfv
    
    #Model Prediction in ln units
    sigma = (tauM**2 + phi**2)**0.5
    
    return sigma
    
    
def bssa14_stdev(M,Rjb,Vs30,intensity_measure='PGA'):
    '''
    Get GMPE standar deviation
    '''
    from numpy import log,exp,sqrt,array,ones,where,zeros
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    
    #Convert input to floats
    Vs30=Vs30.astype('float')
    Rjb=Rjb.astype('float')
    M=M.astype('float')
    
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    elif intensity_measure.upper()=='SA0.3':
        coefficients=[1.2217,1.2401,1.0246,1.2653,0.95676,-0.1959,-0.092855,6.14,-1.09480,0.13388,-0.00548,4.5,
                        1,4.93,0.000000,0.002200,-0.003300,-0.8417,1308.47,760,0,0.1,-0.2191,-0.00670,-9.9,-9.9,
                        103.150,268.590,0.138,0.050,225,300,0.6750,0.5610,0.3630,0.2290]
    elif intensity_measure.upper()=='SA1.0':
        coefficients=[0.3932,0.4218,0.207,0.4124,1.5004,-0.18983,0.17895,6.2,-1.19300,0.10248,-0.00121,4.5,1,5.74,
                        0.000000,0.002920,-0.002090,-1.0500,1109.95,760,0,0.1,-0.1052,-0.00844,0.367,0.208,116.390,
                        270.000,0.098,0.020,225,300,0.5530,0.6250,0.4980,0.2980]
    elif intensity_measure.upper()=='SA3.0':
        coefficients=[-1.1898,-1.142,-1.23,-1.2664,2.1323,-0.04332,0.62694,6.2,-1.21790,0.09764,0.00000,4.5,1,6.93,
                        0.000000,0.002620,-0.001190,-1.0112,922.43,760,0,0.1,-0.0136,-0.00183,1.135,0.516,130.360,
                        195.000,0.088,0.000,225,300,0.5340,0.6190,0.5370,0.3440]
    else:
        print 'ERROR: Unknown intensity measure'
        #return

    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    C = coefficients[17]
    Vc = coefficients[18]
    Vref = coefficients[19]
    f1 = coefficients[20]
    f3 = coefficients[21]
    f4 = coefficients[22]
    f5 = coefficients[23]
    f6 = coefficients[24]
    f7 = coefficients[25]
    R1 = coefficients[26]
    R2 = coefficients[27]
    Dfr = coefficients[28]
    Dfv = coefficients[29]
    V1 = coefficients[30]
    V2 = coefficients[31]
    phi1 = coefficients[32]
    phi2 = coefficients[33]
    tau1 = coefficients[34]
    tau2 = coefficients[35]
    
    tauM=zeros(M.shape)
    #if M <= 4.5:
    i=where(M<=4.5)
    tauM[i] = tau1
    #elif M > 4.5 and M < 5.5:
    i=where((M>4.5) & (M<5.5))[0]
    tauM[i] = tau1 + (tau2 - tau1) * (M[i] - 4.5)
    #else:
    i=where(M>=5.5)[0]
    tauM[i] = tau2
    
    phiM=zeros(M.shape)
    #if M <= 4.5:
    i=where(M<=4.5)
    phiM[i] = phi1
    #elif M > 4.5 and M < 5.5:
    i=where((M>4.5) & (M<5.5))[0]
    phiM[i] = phi1 + (phi2 - phi1) * (M[i] - 4.5)
    #else:
    i=where(M>=5.5)[0]
    phiM[i] = phi2

    phiMR=zeros(M.shape)
    #if Rjb <= R1:
    i=where(Rjb<=R1)[0]
    phiMR[i] = phiM[i]
    #elif Rjb > R1 and Rjb <= R2:
    i=where((Rjb>R1) & (Rjb<=R2))[0]
    phiMR[i] = phiM[i] + Dfr * (log(Rjb[i] / R1) / (log(R2 / R1)))
    #else:
    i=where(Rjb>R2)[0]
    phiMR[i] = phiM[i] + Dfr
    
    phi=zeros(M.shape)
    #if Vs30 >= V2:
    i=where(Vs30>=V2)[0]
    phi[i] = phiMR[i]
    #elif Vs30 >= V1 and Vs30 <= V2:
    i=where((Vs30>=V1) & (Vs30<=V2))[0]
    phi[i] = phiMR[i] - Dfv * (log(V2 / Vs30[i]) / (log(V2 / V1)))
    #else:
    i=where(Vs30<V1)[0]
    phi[i] = phiMR[i] - Dfv
    
    
    #Model Prediction in ln units
    sigma = (tauM**2 + phi**2)**0.5
    
    return sigma
    
def sample_gmpe(median_motion,stdev):
    '''
    Obtain ground motion from lognormal distribution
    '''
    
    from numpy.random import randn
    from numpy import log,exp
    
    x=randn(len(median_motion))
    mean_pga=log(median_motion) #because GMPE's work in Ln space
    #Scale to lognormal
    z_pga=x*stdev+mean_pga
    #And back to physical units
    sampled_motion=exp(z_pga)
    
    return sampled_motion