from numpy import arange,log,exp,r_
from matplotlib import pyplot as plt
from scipy.special import gamma
import Cua2008
from numpy import fft,sin,pi
from numpy.random import normal

duration=60
hf_dt=0.01
mean=0.0
std=1.0
num_samples = int(duration/hf_dt)
t=arange(0,duration,hf_dt)
noise = normal(mean, std, size=num_samples)
freq=0.1
freq2=20.0
noise=sin(2*pi*t*freq+pi/4)+sin(2*pi*t*freq2+pi/6)
ft=fft.rfft(noise)
f=fft.rfftfreq(len(noise),hf_dt)

# GP window
Tw=duration
epsilon=0.2
eta=0.05


b=-epsilon*log(eta)/(1+eta*(log(epsilon)-1))
c=b/(epsilon*Tw)
#a=(exp(1)/(epsilon*Tw))**b
a=(((2*c)**(2*b+1))/gamma(2*b+1))**0.5
w=a*t**b*exp(-c*t)
plt.figure()
plt.plot(r_[0,t+10],r_[0,w])


#fft




#Cua window
i = 0 # Horizontal P-wave acceleration - rock:
i = 6 # Vertical P-wave acceleration - rock:
i = 12 # Horizontal S-wave acceleration - rock: **BAD #<-- S-wave alpha_t_rise for i=12 should be 0.064 instead of 0.64?
i = 19 # Vertical S-wave acceleration - soil:

M=5
R=10
TT=10
env=Cua2008.envelope(M,R,t,TT,Pcoeff=0,Scoeff=12)
plt.figure()
plt.plot(t,env)

plt.show()