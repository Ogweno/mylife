#!/usr/bin/env python

# A module for calculating envelopes based on Cua and Heathon 2008

# By Ran Novitsky Nof @ BSL, 2016

# based on Angie Chung script

#/**********************************************************************************

#*    Copyright (C) by Ran Novitsky Nof                                            *

#*                                                                                 *

#*    This file is part of E2Envelope                                              *

#*                                                                                 *

#*    ElViS is free software: you can redistribute it and/or modify                *

#*    it under the terms of the GNU Lesser General Public License as published by  *

#*    the Free Software Foundation, either version 3 of the License, or            *

#*    (at your option) any later version.                                          *

#*                                                                                 *

#*    This program is distributed in the hope that it will be useful,              *

#*    but WITHOUT ANY WARRANTY; without even the implied warranty of               *

#*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *

#*    GNU Lesser General Public License for more details.                          *

#*                                                                                 *

#*    You should have received a copy of the GNU Lesser General Public License     *

#*    along with this program.  If not, see <http://www.gnu.org/licenses/>.        *

#***********************************************************************************/





from numpy import log10,where,arctan,pi,exp,zeros,sqrt



#Velocity constants (for estimating S wave Trigger time)

Vp,Vs = 5.6,3.5



######################################################################

        # Envelope equation from Cua thesis, 2004; Cua and Heaton, 2008 (submitted to BSSA):



        # Coefficients:

        #i = 0 # Horizontal P-wave acceleration - rock:

        #i = 1 # Horizontal P-wave acceleration - soil:  **BAD #<-- P-wave alpha_gamma for i=1 should be =-0.048, not =-0.48?

        #i = 2 # Horizontal P-wave velocity - rock:

        #i = 3 # Horizontal P-wave velocity - soil:

        #i = 4 # Horizontal P-wave (filtered) displacement - rock:

        #i = 5 # Horizontal P-wave (filtered) displacement - soil:

        #i = 6 # Vertical P-wave acceleration - rock:

        #i = 7 # Vertical P-wave acceleration - soil:

        #i = 8 # Vertical P-wave velocity - rock:

        #i = 9 # Vertical P-wave velocity - soil:

        #i = 10 # Vertical P-wave (filtered) displacement - rock:

        #i = 11 # Vertical P-wave (filtered) displacement - soil:

        #i = 12 # Horizontal S-wave acceleration - rock: **BAD #<-- S-wave alpha_t_rise for i=12 should be 0.064 instead of 0.64?

        #i = 13 # Horizontal S-wave acceleration - soil:

        #i = 14 # Horizontal S-wave velocity - rock:

        #i = 15 # Horizontal S-wave velocity - soil:

        #i = 16 # Horizontal S-wave (filtered) displacement - rock:

        #i = 17 # Horizontal S-wave (filtered) displacement - soil:

        #i = 18 # Vertical S-wave acceleration - rock:

        #i = 19 # Vertical S-wave acceleration - soil:

        #i = 20 # Vertical S-wave velocity - rock:

        #i = 21 # Vertical S-wave velocity - soil:

        #i = 22 # Vertical S-wave (filtered) displacement - rock:

        #i = 23 # Vertical S-wave (filtered) displacement - soil:



Pcoeff = 6 # Vertical P-wave acceleration - rock:

Scoeff = 18 # Vertical S-wave acceleration - rock:



a = [0.719, 0.737, 0.801, 0.836, 0.950, 0.943, 0.745, 0.739, 0.821, 0.812, 0.956, 0.933,

            0.779, 0.836, 0.894, 0.960, 1.031, 1.081, 0.778, 0.751, 0.900, 0.882, 1.042, 1.034]

b = [-3.273e-3, -2.520e-3, -8.397e-4, -5.409e-4, -1.685e-6, -5.171e-7, -4.010e-3, -4.134e-3,

            -8.543e-4, -2.652e-6, -1.975e-6, -1.090e-7, -2.555e-3, -2.324e-3, -4.286e-4, -8.328e-4,

            -1.015e-7, -1.204e-6, -2.66e-5, -2.473e-3, -1.027e-5,- 5.41e-4, -1.124e-5, -4.924e-6]

d = [-1.195, -1.26, -1.249, -1.284, -1.275, -1.161, -1.200, -1.199, -1.362, -1.483, -1.345, -1.234,

            -1.352, -1.562, -1.440, -1.589, -1.438, -1.556, -1.385, -1.474, -1.505, -1.484, -1.367, -1.363]

c1 = [1.600, 2.410, 0.761, 1.214, 2.162, 2.266, 1.752, 2.030, 1.148, 1.402, 1.656, 1.515,

            1.478, 2.423, 1.114, 1.982, 1.098, 1.946, 1.763, 1.593, 1.388, 1.530, 1.379, 1.549]

c2 = [1.045, 0.955, 1.340, 0.978, 1.088, 1.016, 1.091, 1.972, 1.100, 0.995, 1.164, 1.041,

            1.105, 1.054, 1.110, 1.067, 1.133, 1.091, 1.112, 1.106, 1.096, 1.04, 1.178, 1.082]

e = [-1.065, -1.051, -3.103, -3.135, -4.958, -5.008, -0.955, -0.775, -2.901, -2.551, -4.799, -4.749,

            -0.645, -0.338, -2.602, -2.351, -4.342, -4.101, -0.751, -0.355, -2.778, -2.537, -4.738, -4.569]

sig_uncorr = [0.307, 0.286, 0.268, 0.263, 0.284, 0.301, 0.288, 0.317, 0.263, 0.298, 02.83, 0.312,

            0.308, 0.312, 0.279, 0.296, 0.277, 0.326, 0.300, 0.300, 0.250, 0.270, 0.253, 0.286]

sig_corr = [0.233, 0.229, 0.211, 0.219, 0.239, 0.247, 0.243, 0.256, 0.231, 0.239, 0.254, 0.248,

            0.243, 0.248, 0.230, 0.230, 0.233, 0.236, 0.238, 0.235, 0.220, 0.221, 0.232, 0.230]



        # Coefficients for eqn: log(env_param) = alpha*M + beta*R + delta*logR + mu

        # Coefficients and equation for t_rise (rise time):

alpha_t_rise = [0.06, 0.07, 0.06, 0.07, 0.05, 0.05, 0.06, 0.06, 0.06, 0.06, 0.08, 0.067,

            0.064, 0.055, 0.093, 0.087, 0.109, 0.12, 0.069, 0.059, 0.116, 0.11, 0.123, 0.124]  #<- typo in i=12? 0.064 instead of 0.64?

beta_t_rise = [5.5e-4, 1.2e-3, 1.33e-3, 4.35e-4, 1.29e-3, 1.19e-3, 7.45e-4, 5.87e-4, 7.32e-4, 1.08e-3, 1.64e-3, 1.21e-3,

            0, 1.21e-3, 0, 4.0e-4, 7.68e-4, 0, 0, 2.18e-3, 0, 1.24e-3, 1.3e-3, 0]

delta_t_rise = [0.27, 0.24, 0.23, 0.47, 0.27, 0.47, 0.37, 0.23, 0.25, 0.22, 0.13, 0.28,

            0.48, 0.34, 0.48, 0.49, 0.38, 0.45, 0.49, 0.26, 0.503, 0.38, 0.257, 0.439]

mu_t_rise = [-0.37, -0.38, -0.34, -0.68, -0.34, -0.58, -0.51, -0.37, -0.37, -0.36, -0.33, -0.46,

            -0.89, -0.66, -0.96, -0.98, -0.87,-0.89,-0.97, -0.66, -1.14, -0.91, -0.749, -0.82]



        # Coefficients and equation for delta_t (wave duration):

alpha_delta_t = [0, 0.03, 0.054, 0.03, 0.047, 0.051, 0, 0, 0.046, 0.031, 0.058, 0.043,

            0, 0.028, 0.02, 0.028, 0.04, 0.03, 0.03, 0.03, 0.018, 0.017, 0.033, 0.023]

beta_delta_t = [2.58e-3, 2.37e-3, 1.93e-3, 2.03e-3, 0, 1.12e-3, 2.75e-3, 1.76e-3, 2.61e-3, 1.7e-3, 2.02e-3, 9.94e-4,

            -4.87e-4, 0, 0, 0, 1.1e-3, 0, -1.4e-3, -1.78e-3, 0, -6.93e-4, 2.6e-4, -7.18e-4]

delta_delta_t = [0.21, 0.39, 0.16, 0.289, 0.45, 0.33, 0.165, 0.36, 0, 0.26, 0, 0.19,

            0.13, 0.07, 0, 0.046, -0.15, 0.037, 0.22, 0.307, 0, 0.119, 0, 0.074]

mu_delta_t = [-0.22, -0.59, -0.36, -0.45, -0.68, -0.59, -0.245, -0.48, -0.213, -0.52, -0.253, -0.42,

            0.0024, -0.102, 0.046, -0.083, 0.11, -0.066, -0.17, -0.66, -0.072, -0.05, -0.015, -0.005]



        # Coefficients and equation for tau (decay):

alpha_tau = [0.047, 0.087, 0.054, 0.0403, 0, 0.035, 0.03, 0.057, 0.03, 0.0311, 0.05, 0.052,

            0.037, 0.0557, 0.029, 0.045, 0.029, 0.038, 0.031, 0.06, 0.04, 0.051, 0.024, 0.022]  #<- typo in i=9? 0.031 instead of 0.31?

beta_tau = [0, -1.89e-3, 5.37e-5, -1.26e-3, 0, -1.27e-3, 2.75e-3, -1.36e-3, 8.6e-4, -6.4e-4, 8.9e-4, 0,

            0, -8.2e-4, 8.0e-4, -5.46e-4, 0, -1.34e-3, 0, -1.45e-3, 9.4e-4, -1.41e-3, 0, -1.65e-3]

delta_tau = [0.48, 0.58, 0.41, 0.387, 0.19, 0.19, 0.58, 0.63, 0.35, 0.44, 0.16, 0.12,

            0.39, 0.51, 0.25, 0.46, 0.36, 0.48, 0.34, 0.51, 0.25, 0.438, 0.303, 0.44]

gamma_tau = [0.82, 0.58, 0.73, 0.58, 0, 0, 0, 0, 0, 0, 0, 0, 1.73, 1.63, 1.61, 0, 0, 0, 0, 0, 0, 0, 0, 0]

mu_tau = [-0.75, -0.87, -0.51, -0.372, -0.07, -0.03, -0.97, -0.96, -0.62, -0.55, -0.387, -0.166,

            -0.59, -0.68, -0.31, -0.55, -0.38, -0.39, -0.44, -0.60, -0.34, -0.368, -0.22, -0.19]

avg_gamma = 0.15



# Coefficients and equation for gamma (decay):

alpha_gamma = [-0.032, -0.048, -0.044, -0.0403, -0.062, -0.061, -0.027, -0.024, -0.039, -0.037, -0.052, -0.066,

            -0.014, -0.015, -0.024, -0.031, -0.025, -2.67e-2, -0.0149, -0.0197, -0.028, -0.0334, -0.015, -0.0176] #<--should be =-0.048 for i=1? not =-0.48?

beta_gamma = [-1.81e-3, -1.42e-3, -1.65e-3, -2.0e-3, -2.3e-3, -1.9e-3, -1.75e-3, -1.6e-3, -1.88e-3, -2.23e-3, -1.67e-3, -2.5e-3,

            -5.28e-4, -5.89e-4, -1.02e-3, -4.61e-4, -4.22e-4, 2.0e-4, -4.64e-4, 0, -8.32e-4, 0, 0, 5.65e-4]

delta_gamma = [-0.1, -0.13, -0.16, 0, 0, 0.11, -0.18, -0.24, -0.18, -0.14, -0.21, 0,

            -0.11, -0.163, -0.055, -0.162, -0.145, -0.217, -0.122, -0.242, -0.123, -0.21, -0.229, -0.25]

tau_gamma = [0.27, 0.26, 0.33, 0, 0, 0, 0, 0, 0, 0, 0, 0,

            0.38, 0.39, 0.36, 0, 0, 0, 0, 0, 0, 0, 0, 0]

mu_gamma = [0.64, 0.71, 0.72, 0.578, 0.61, 0.39, 0.74, 0.84, 0.76, 0.71, 0.849, 0.63,

            0.26, 0.299, 0.207, 0.302, 0.262, 0.274, 0.255, 0.378, 0.325, 0.325, 0.309, 0.236]

avg_gamma = 0.15



# place holders for error estimators and station corrections

stat_err = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

sta_corr =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]



def envelope(M,R,times,TT,Pcoeff=Pcoeff,Scoeff=Scoeff):

  '''

  Calculates Envelopes.

  INPUTS:

    M - magnitude as float

    R - Distance in km

    times - time series as floats (in seconds)

    TT - trigger time as float (in seconds)

    Pcoeff, Scoeff - coefficient index for P and S waves (see table at top). default to vertical acceleration - rock

  '''

  TTS = TT+R/Vs # S arrival


  # coefficients

  t_rise_p = 10**(alpha_t_rise[Pcoeff] * M + beta_t_rise[Pcoeff] * R + delta_t_rise[Pcoeff] * log10(R) + mu_t_rise[Pcoeff])

  t_rise_s = 10**(alpha_t_rise[Scoeff] * M + beta_t_rise[Scoeff] * R + delta_t_rise[Scoeff] * log10(R) + mu_t_rise[Scoeff])

  delta_t_p = 10**(alpha_delta_t[Pcoeff] * M + beta_delta_t[Pcoeff] * R + delta_delta_t[Pcoeff] * log10(R) + mu_delta_t[Pcoeff])

  delta_t_s = 10**(alpha_delta_t[Scoeff] * M + beta_delta_t[Scoeff] * R + delta_delta_t[Scoeff] * log10(R) + mu_delta_t[Scoeff])

  tau_p = 10**(alpha_tau[Pcoeff] * M + beta_tau[Pcoeff] * R + delta_tau[Pcoeff] * log10(R) + mu_tau[Pcoeff])

  tau_s = 10**(alpha_tau[Scoeff] * M + beta_tau[Scoeff] * R + delta_tau[Scoeff] * log10(R) + mu_tau[Scoeff])

  gamma_p = 10**(alpha_gamma[Pcoeff] * M + beta_gamma[Pcoeff] * R + delta_gamma[Pcoeff] * log10(R) + mu_gamma[Pcoeff])

  gamma_s = 10**(alpha_gamma[Scoeff] * M + beta_gamma[Scoeff] * R + delta_gamma[Scoeff] * log10(R) + mu_gamma[Scoeff])

  # Other variable (turn on saturation for larger events?)

  C_p = (arctan(M-5) + (pi/2))*(c1[Pcoeff]*exp(c2[Pcoeff] * (M-5)))

  C_s = (arctan(M-5) + (pi/2))*(c1[Scoeff]*exp(c2[Scoeff] * (M-5)))

  R1 = sqrt(R**2 + 9)

  # Basic AMplitudes

  A_p = 10**(a[Pcoeff]*M + b[Pcoeff]*(R1 + C_p) + d[Pcoeff]*log10(R1+C_p) + e[Pcoeff]+(sta_corr[Pcoeff]) + stat_err[Pcoeff])

  A_s = 10**(a[Scoeff]*M + b[Scoeff]*(R1 + C_s) + d[Scoeff]*log10(R1+C_s) + e[Scoeff]+(sta_corr[Scoeff]) + stat_err[Scoeff])

  # calculate envelope (ENV)

  ENV = zeros(len(times))

  # P envelope

  indx = where((times>=TT) & (times<TT+t_rise_p)) # between trigger and rise time

  if len(indx): ENV[indx] = (A_p/t_rise_p*(times[indx]-TT)) # make sure we have data in that time frame and get envelope

  indx = where((times>=TT+t_rise_p) & (times<TT+t_rise_p+delta_t_p)) # flat area

  if len(indx): ENV[indx] = A_p # make sure we have data in that time frame and get envelope

  indx = where(times>TT+t_rise_p+delta_t_p) # coda

  if len(indx): ENV[indx] = (A_p/((times[indx]-TT-t_rise_p-delta_t_p+tau_p)**gamma_p)) # make sure we have data in that time frame and get envelope

  # S envelope

  indx = where((times>=TTS) & (times<TTS+t_rise_s)) # between trigger and rise time

  if len(indx): ENV[indx] += (A_s/t_rise_s*(times[indx]-TTS)) # make sure we have data in that time frame and get envelope

  indx = where((times>=TTS+t_rise_s) & (times<TTS+t_rise_s+delta_t_s)) # flat area

  if len(indx): ENV[indx] += A_s # make sure we have data in that time frame and get envelope

  indx = where(times>TTS+t_rise_s+delta_t_s) # coda

  if len(indx): ENV[indx] += (A_s/((times[indx]-TTS-t_rise_s-delta_t_s+tau_s)**gamma_s)) # make sure we have data in that time frame and get envelope

  return ENV # envelope data of P and S calculated at the same time points as input times variable.







