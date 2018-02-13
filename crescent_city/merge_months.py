from obspy import read

#st=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_01_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_02_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_03_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_04_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_05_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_06_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_07_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_08_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_09_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_10_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_11_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_12_2017.sac')
#
#st.merge()
#st.write(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_2017.sac',format='SAC')


#st=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_01_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_02_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_03_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_04_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_05_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_06_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_07_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_08_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_09_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_10_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_11_2017.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/eureka/eure_12_2017.sac')
#
#st.merge(method=1)
#st.write(u'/Users/dmelgar/tidegauge_noise/eureka/eure_2017.sac',format='SAC')


#st=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_10_2012.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_11_2012.sac')
#
#
#st.merge(method=1)
#st.write(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_1011_2012.sac',format='SAC')





#st=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_09_2009.sac')
#st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_10_2009.sac')
#
#
#st.merge(method=1)
#st.write(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_0910_2009.sac',format='SAC')




st=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_02_2010.sac')
st+=read(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_03_2010.sac')


st.merge(method=1)
st.write(u'/Users/dmelgar/tidegauge_noise/crescent_city/cres_0203_2010.sac',format='SAC')