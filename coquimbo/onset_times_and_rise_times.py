from numpy import genfromtxt,deg2rad,cos,sin,savetxt
from mudpy import fakequakes as fq

f=genfromtxt(u'/Users/dmelgar/Coquimbo2015/coquimbo_FQ_3km.rupt')
fmoment=genfromtxt('/Volumes/Illapel/FQ/illapel/output/ruptures/illapel.000000.rupt',usecols=13)
fout=u'/Users/dmelgar/Coquimbo2015/coquimbo_FQ_3km_stoc.0000.rupt'
home='/Volumes/Illapel/FQ/'
project_name='illapel'
run_name='illapel'
model_name='maule.mod' 
hypocenter=[-71.654,-31.570,29.8]

#Get stochastic rake vector
stoc_rake=fq.get_stochastic_rake(90,len(f))

slip=f[:,9]
Mw=8.3
M0=10**(1.5*Mw+9.1)
rise_times=fq.get_rise_times(M0,slip,f,[5,8],stoc_rake)
t_onset=fq.get_rupture_onset(home,project_name,slip,f,model_name,hypocenter,[5,8],M0)

f[:,8]=slip*cos(deg2rad(stoc_rake))
f[:,9]=slip*sin(deg2rad(stoc_rake))
f[:,7]=rise_times
f[:,12]=t_onset
f[:,13]=fmoment

#   1	-72.5799	-30.0132	7.8040	3.19	6.65	0.50	5.00	0.00	5.59	2692.58	2692.58	0.00	0.00
savetxt(fout,f,fmt='%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t')
