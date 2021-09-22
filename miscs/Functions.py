#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 22:47:46 2019

@author: giordano
"""

import numpy as np
import math
from scipy import signal
import scipy.fftpack as fftpack
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt



#------------------------------------------------------------------------------
# F U N C T I O N   F O R   L I N E A R   R E G R E S S I O N   C O E F F
#------------------------------------------------------------------------------

# N - number of days after which surface runoff ceases, A should be in km**2
def N_g(A):
    # N - number of days after which surface runoff ceases, A should be in km**2
    N = .83*A**(.2)  
    # N2 - the odd integer between 3 and 11 nearest to 2N    
    if int(2*N)%2 == 0:
        N2 = int(2*N) + 1
    else:
        N2 = int(2*N)       
    return N2

def func_regr(x, alpha):
    return alpha*x

# Define a routin to calculate alpha
def alpha_reg(data,A):
    # Number of days for regression
    Ng = N_g(A)   
    # i_decr is an array for the index of decreasing discharges
    i_decr = np.array([])
    for i in range(1,len(data)):
        if data[i] < data[i-1]:
            i_decr = np.append(i_decr,i)
    # transform the index in int            
    i_decr = i_decr.astype(int)
    
    # i_sep contains the index of the fist and the last element of decreasing discharge set
    i_sep = np.array([i_decr[0]])
    for i in range(1,len(i_decr)):
        if i_decr[i] != i_decr[i-1] + 1:
            i_sep = np.append(i_sep,np.array([i_decr[i-1],i_decr[i]]))
    # add to i_sep the last element of the entire series
    i_sep = np.append(i_sep,i_decr[-1])
    # reshape the array such that [[start1, end1],[start2, end2]...]
    i_sep = np.reshape(i_sep,(int(len(i_sep)/2),2))
    
    # index_rem is an array containing the index of the set having less than 5 elements and so needed be deleted
    index_rem = np.array([])    
    for i in range(len(i_sep)):
        if i_sep[i,1] < i_sep[i,0] + Ng:
            index_rem = np.append(index_rem,i)          
    index_rem = index_rem.astype(int)
    
    # detete from i_sep the index of the decreasing index having less than 5 elements
    i_sep_star = np.delete(i_sep, index_rem, 0)
    
    # creation of two arrays: discharge for time = tau and discharge for time = tau + 1    
    N_delay = int(N_g(A)/4) # Days after which calculating the recession curve
    btau = np.array([])
    btau1 = np.array([])
    for i in range(len(i_sep_star)):
        btau = np.append(btau,data[i_sep_star[i,0]+N_delay:i_sep_star[i,1]])
        btau1 = np.append(btau1,data[i_sep_star[i,0]+N_delay+1:i_sep_star[i,1]+1])
    
    # parameters for linear regression
    alpha, pcov = curve_fit(func_regr, btau, btau1, bounds=(0, np.inf), maxfev=1000)
    return alpha

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# B A S E F L O W   F O R   E C K H A R D T
#------------------------------------------------------------------------------

# Define the filter for the Eckman algorithm
def filter_Eck(BFImax, data, A):
    # Calculate alpha
    alpha = alpha_reg(data, A)
    res = np.zeros(len(data))
    res[0] = ((1-alpha)*BFImax*data[0])/(1-alpha*BFImax)
    for i in range(1,len(data)):
        res[i] = ((1-BFImax)*alpha*res[i-1]+(1-alpha)*BFImax*data[i])/(1-alpha*BFImax)
    # If the baseflow is larger than the measured discharge, then it is equal to the discharge
    for i in range(len(data)):
        if res[i] > data[i]: res[i] = data[i]
    return res


def BF_Eck(data, A):
    # Define the parameter BFImax
    BFImax = 0.80    
    
    # First pass 
    QB1 = filter_Eck(BFImax, data, A)
    # Second pass
    QB2 = np.flip(filter_Eck(BFImax, np.flip(QB1, axis=None), A))
    # Third pass
    QB = filter_Eck(BFImax, QB2, A)

    return QB


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# B A S E F L O W   F O R   A R N O L D
#------------------------------------------------------------------------------

def filter_LyneHollick(alpha, data):
    # Fist pass filter
    res_fast = np.zeros(len(data))
    res = np.zeros(len(data))
    
    res_fast[0] = (data[0])/2.
    res[0] = (data[0])/2.
    
    for i in range(1,len(data)):
        res_fast[i] = (alpha*res_fast[i-1] + (1+alpha)/2.*(data[i]-data[i-1]))
        res[i] = data[i] - res_fast[i]
        if res[i] > data[i]: res[i] = data[i] 
    return res

def BF_LyneHollick(data):

    # Filter parameter
    alpha = 0.925
    
    # First pass 
    QB1 = filter_LyneHollick(alpha, data)
    # Second pass
    QB2 = np.flip(filter_LyneHollick(alpha, np.flip(QB1, axis=None)))
    # Third pass
    QB = filter_LyneHollick(alpha, QB2)
    return QB

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# B A S E F L O W   F O R   C H A P M A N
#------------------------------------------------------------------------------

def BF_Chapman(data, A):
    
    alpha = alpha_reg(data, A)
    
    # Fist pass filter
    qt = np.zeros(len(data))
    QB = np.zeros(len(data))
    
    qt[0] = (data[0])/2.
    QB[0] = (data[0])/2.
    
    
    for i in range(1,len(data)):
        qt[i] = ((3*alpha-1)/(3-alpha)*qt[i-1] + (2)/(3-alpha)*(data[i]-alpha*data[i-1]))
        QB[i] = data[i] - qt[i]
        if QB[i] > data[i]: QB[i] = data[i]
    
    return QB

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# B A S E F L O W   F O R   U K H I
#------------------------------------------------------------------------------

def BF_UKHI(data, A):
    
    # Set the number of days within calculate the minumum
    Ng = N_g(A)   
    
    # Set the array containing the index of the minimum in Ng days
    Qimin = np.array([])
    # Find the index for the local minimum in each i-th interval of 5-days discharge    
    for i in range(math.ceil(len(data)/Ng)):
        start = i*Ng
        end = i*Ng+(Ng-1)
        if end <= len(data):
            Q5 = data[start : end+1]
        else:
            Q5 = data[start : len(data)]
        Qimin = np.append(Qimin, np.where(Q5 == min(Q5))[0] + start)
    # Convert the float array in a int array        
    Qimin = Qimin.astype(int)
    
    # Initialized the vector containing the index for the turning point (baseflow points) following the rule 0.9*Q[i] < [Q[i-1],Q[i+1]]
    QBindex = np.array([])  
    # Find the index of the turning points  
    for i in range(1,len(Qimin)-1):
        if (0.9*data[Qimin[i]] <= data[Qimin[i-1]] and 0.9*data[Qimin[i]] <= data[Qimin[i+1]]):
            QBindex = np.append(QBindex, Qimin[i])
    # Convert the QB index from float to integer        
    QBindex = QBindex.astype(int)
 
    # Initialize the vector for the Base Flow
    QB = np.zeros(len(data))

    # The baseflow is equal to the discharge for the turning points
    QB[QBindex] = data[QBindex]

    # Define a vector containing the angular coefficient between local minima
    linearizSector = np.zeros(len(QBindex) + 1)

    # Calculate the local minimum    
    for i in range(1,len(linearizSector)-1):
        linearizSector[i] =(data[QBindex[i]]-data[QBindex[i-1]])/(QBindex[i] - QBindex[i-1])

    # Calculate the baseflow multiplying the angular coefficient with the position    
    for i in range(1,len(linearizSector)-1):
        for j in range(QBindex[i-1]+1,QBindex[i]):
            QB[j] = QB[QBindex[i-1]] + linearizSector[i]*(j-QBindex[i-1])

    # If the the local minimum is not at zero, calculate the BF for the previous indeces        
    if QBindex[0] != 0:
        for i in range(QBindex[0]):
            QB[i] = i*data[QBindex[0]]/QBindex[0]

    # If the the local minimum is not at the end of observations, calculate the BF for the sequent indeces            
    if QBindex[-1] != (len(data) - 1):
        for i in range(QBindex[-1]+1,len(data)):
            QB[i] = data[QBindex[-1]] + data[QBindex[-1]]/(len(data)-QBindex[-1])*(QBindex[-1]-i)
    
    # If the linearized BF if larger than the measured discharge, than BF = measured discharge            
    for i in range(len(QB)):
        if QB[i] > data[i]: QB[i] = data[i]
            
    return QB

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# B A S E F L O W   F O R   H Y S E P
#------------------------------------------------------------------------------
# HYSEP 1 - FIXED-INTERVAL METHOD

def BF_HYSEP1(data, A):
    # N - number of days after which surface runoff ceases, A should be in km**2
    N = .83*A**(.2)  

    # N2 - the odd integer between 3 and 11 nearest to 2N    
    if int(2*N)%2 == 0:
        N2 = int(2*N) + 1
    else:
        N2 = int(2*N)
    
    # Set to zero the QB vector with the same length of baseflow    
    QB = np.zeros(len(data))
    
    # find the lowest discharge in each interval of N2 and assing it to the baseflow in the same interval
    for i in range(0,int(len(QB)/N2)*N2,N2):
        QB[i : i + N2] = min(data[i : i + N2])
    
    if int(len(QB)/N2)*N2 < len(QB):
        QB[int(len(QB)/N2)*N2:] = min(data[int(len(QB)/N2)*N2:])   
    
    return QB
    
#------------------------------------------------------------------------------
# HYSEP 2 - SLIDING-INTERVAL METHOD

def BF_HYSEP2(data, A):
    # N - number of days after which surface runoff ceases, A should be in km**2
    N = .83*A**(.2)  

    # N2 - the odd integer between 3 and 11 nearest to 2N        
    if int(2*N)%2 == 0:
        N2 = int(2*N) + 1
    else:
        N2 = int(2*N)
    
    # Set to zero the QB vector with the same length of baseflow        
    QB = np.zeros(len(data))

    # The sliding-interval method finds the lowest discharge in one half the interval minus 1 day [0.5(2N*-1) days] before and after the day being considered and assigns it to that day    
    for i in range(0,len(QB)):
        ileft = i-int(0.5*(N2-1))
        if ileft < 0: ileft = 0
        iright = i+int(0.5*(N2-1))
        if iright > len(QB) : iright = len(QB)
        QB[i] = min(data[ileft:iright])

    return QB

#------------------------------------------------------------------------------
# HYSEP 3 - LOCAL-MINIMUM METHOD

def BF_HYSEP3(data, A):
    # N - number of days after which surface runoff ceases, A should be in km**2
    N = .83*A**(.2)  

    # N2 - the odd integer between 3 and 11 nearest to 2N            
    if int(2*N)%2 == 0:
        N2 = int(2*N) + 1
    else:
        N2 = int(2*N)
    
    # Set the array collecting the index of local minimum    
    QBindex = np.array([])
  
    #The local-minimum method checks each day to determine if it is the lowest discharge in one half the interval minus 1 day [0.5(2N*-1) days] before and after the day being considered.
    for i in range(0,len(data)):
        # ileft and iright are the indeces of the interval centered in i
        ileft = i-int(0.5*(N2-1))
        if ileft < 0: ileft = 0
        iright = i+int(0.5*(N2-1))
        if iright > len(data) : iright = len(data)
        
        # if the discharge in the i-th instant is equal to the minimim between ileft and i right, the discharge in i is baseflow
        if data[i] == min(data[ileft:iright]):
            QBindex = np.append(QBindex, i)
    
    # QBindex is a vector of integer    
    QBindex = QBindex.astype(int)
    
    # Set QB as a vector of zero of the same length of discharge
    QB = np.zeros(len(data))

    # QB in the QBindex is equal to the discharge in QBindex
    QB[QBindex] = data[QBindex]

    # Define a vector containing the angular coefficient between local minima
    linearizSector = np.zeros(len(QBindex) + 1)

    # Calculate the local minimum    
    for i in range(1,len(linearizSector)-1):
        linearizSector[i] =(data[QBindex[i]]-data[QBindex[i-1]])/(QBindex[i] - QBindex[i-1])

    # Calculate the baseflow multiplying the angular coefficient with the position    
    for i in range(1,len(linearizSector)-1):
        for j in range(QBindex[i-1]+1,QBindex[i]):
            QB[j] = QB[QBindex[i-1]] + linearizSector[i]*(j-QBindex[i-1])

    # If the the local minimum is not at zero, calculate the BF for the previous indeces        
    if QBindex[0] != 0:
        for i in range(QBindex[0]):
            QB[i] = i*data[QBindex[0]]/QBindex[0]

    # If the the local minimum is not at the end of observations, calculate the BF for the sequent indeces            
    if QBindex[-1] != (len(data) - 1):
        for i in range(QBindex[-1]+1,len(data)):
            QB[i] = data[QBindex[-1]] + data[QBindex[-1]]/(len(data)-QBindex[-1])*(QBindex[-1]-i)
    
    # If the linearized BF if larger than the measured discharge, than BF = measured discharge            
    for i in range(len(QB)):
        if QB[i] > data[i]: QB[i] = data[i]

        
    return QB

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Define a function that calculates the BF according the the key word describing the BF technique
    
def BaseFlow(BF_filter, data, A):
    if BF_filter == 'Eckhardt':
        BF = BF_Eck(data, A)
    if BF_filter == 'LyneHollick':
        BF = BF_LyneHollick(data)
    if BF_filter == 'Chapman':
        BF = BF_Chapman(data, A)
    if BF_filter == 'UKHI':
        BF = BF_UKHI(data, A)
    if BF_filter == 'HYSEP1':
        BF = BF_HYSEP1(data, A)
    if BF_filter == 'HYSEP2':
        BF = BF_HYSEP2(data, A)
    if BF_filter == 'HYSEP3':
        BF = BF_HYSEP3(data, A)

    return BF


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# C A L C U L A T I O N   O F   B F I
#------------------------------------------------------------------------------


def BFI(time, discharge, baseflow):
    return np.trapz(baseflow, x=time)/np.trapz(discharge, x=time)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# S P E C T R U M   A N D   P A R A M E T E R S
#------------------------------------------------------------------------------


# function for the spectrum of a time series
def Spectrum_series_fft(fr,vec_val):
    # fr is the frequency of data, vec_val is the array of which we want to calculate spectrum
    spectrum = fftpack.fft(vec_val)
    spectrum = abs(spectrum[: int(round(len(spectrum) / 2))]) ** 2
    spectrum = spectrum[1:]
    frequency = (abs(fftpack.fftfreq(len(vec_val), 1./fr))[: int(round(len(vec_val) / 2))])[1:]    
    # Convert the frequency in angular frequency
    frequency = frequency*2*np.pi
    return (frequency, spectrum)

# function for the spectrum of a time series
def Spectrum_series(fr,vec_val): 
    # fr is the frequency of data, vec_val is the array of which we want to calculate spectrum
    nperseg = int(round(len(vec_val) / 10))
    f_vec, Pxx_vec = signal.welch(vec_val, fr, nperseg=nperseg, window="hamming")
    # convert frequency in angular frequency
    w_vec = f_vec*2*np.pi
    return (w_vec, Pxx_vec)

# Function to obtain a regular array for the spectrum in the log-scale    
#------------------------------------------------------------------------------

def log_scale_spectrum(w_QB, Pxx_QB):
    w_QB_mod = w_QB
#    w_QB_mod[0] = 10**(np.log10(w_QB_mod)[1]-1)

    Pxx_QB_func = interp1d(np.log10(w_QB_mod[1:]), np.log10(Pxx_QB[1:]))
    w_QB_int = np.linspace(min(np.log10(w_QB_mod[1:])), max(np.log10(w_QB_mod[1:])), num=len(Pxx_QB), endpoint=True)
    Pxx_QB_int = Pxx_QB_func(w_QB_int)
    return (w_QB_int, Pxx_QB_int)
    

#------------------------------------------------------------------------------

# Define the coefficient of determination to quantify the goodness of fit
def Rsquared(obs,fit):
    # residual sum of squares
    ss_res = np.sum((obs - fit) ** 2)

    # total sum of squares
    ss_tot = np.sum((obs - np.mean(obs)) ** 2)

    # r-squared
    return 1. - (ss_res / ss_tot)

#------------------------------------------------------------------------------
# Calculation of t_c from variance
#------------------------------------------------------------------------------

def Area_series(X_vec, Y_vec):
    sum_vec = ((X_vec[0]-X_vec[1])/2.)*Y_vec[0]
    for i in range(1,len(X_vec)-1):
        sum_vec = sum_vec + ((X_vec[i]-X_vec[i-1])/2.+(X_vec[i+1]-X_vec[i])/2.)*Y_vec[i]
    sum_vec = sum_vec + ((X_vec[-1]-X_vec[-2])/2.)*Y_vec[-1]
    return sum_vec

# Function to calculate the integral of a spectrum (coinciding with variance)
def Integral_Spectrum(w_QB, Pxx_QB):
    res = Pxx_QB/(2*np.pi)
    return Area_series(w_QB, res)

# Function to calculate the variance of baseflow from the spectrum of recharge
    
def var_qq(Srr,w_rr,t_c):
    t_c_w = t_c/(2.*np.pi)
    Sqq = [Srr[0]/(1+t_c_w**2*w_rr[0]**2)/(2*np.pi)]
    for i in range(1,len(w_rr)):
        Sqq.append(Srr[i]/(1+t_c_w**2*w_rr[i]**2)/(2*np.pi))
    return Area_series(w_rr, Sqq)

# Function to calculate the t_c from the spectrum of Recharge and BF using the variance
# Srr is the spectrum of Recharge, QB is the baseflow timeseries and w_rr is the frequency
def var_VS_tc(Srr,w_rr,tc_max):
    sigmaB = []
    t_c_vec = []
    for t_c in range(tc_max):
        t_c_vec.append(t_c)
        sigmaB.append(var_qq(Srr,w_rr,t_c))
    return (t_c_vec,sigmaB)
    
#------------------------------------------------------------------------------
# Analytical spectrum with two line
#------------------------------------------------------------------------------

def func_tail(x, m, q):
    return m*x + q

def func_smallf(x, m, q):
    return m*x + q

def line_fit_spectrum(w_QB, Pxx_QB):
    (w_QB_int, Pxx_QB_int) = log_scale_spectrum(w_QB, Pxx_QB)
    max_index = max(np.where(Pxx_QB_int>Pxx_QB_int[0]-.5)[0])
    popt_smallf, pcov_smallf = curve_fit(func_smallf, w_QB_int[:max_index], Pxx_QB_int[:max_index])
    popt_tail, pcov_tail = curve_fit(func_tail, w_QB_int[max_index:], Pxx_QB_int[max_index:])
    
    breaking_point = (popt_smallf[1]-popt_tail[1])/(popt_tail[0]-popt_smallf[0])
    breaking_point_tail = np.where(w_QB_int>=breaking_point)
    breaking_point_smallf = np.where(w_QB_int<=breaking_point)
    t_c = 2*np.pi/10**(breaking_point)

    line_spectrum = np.append(func_smallf(w_QB_int[breaking_point_smallf], *popt_smallf),func_tail(w_QB_int[breaking_point_tail], *popt_tail))
    R = Rsquared(Pxx_QB_int,line_spectrum)
    
    return {'analytical_spectrum': np.array([w_QB_int,line_spectrum]), 'Parameters':{'t_c':t_c, 'm1':popt_smallf[0], 'm2':popt_tail[0], 'R':R}}


#------------------------------------------------------------------------------
# Analytical spectrum for linear reservoir
#------------------------------------------------------------------------------

# function of the analytical spectrum as a function of precipitation recharge
def func_wn(w, t_c):
    return 1./(1+(t_c/(2.*np.pi))**2*w**2)

# function of the analytical spectrum as a function of precipitation recharge in log scale
def func_wn_log(w, t_c):
    return np.log10(func_wn(10**w, t_c))

def analytical_spectrum_wn(fr, w_QB, Pxx_QB):
    # spectrum baseflow
    (w_QB_log, Pxx_QB_log) = log_scale_spectrum(w_QB, Pxx_QB)
    
    popt, pcov = curve_fit(func_wn_log, w_QB_log, Pxx_QB_log, bounds=(0, np.inf), maxfev=1000)
    # Paramters fitted
    t_c = popt[0]

    return {'t_c':t_c}

#------------------------------------------------------------------------------
# Analytical spectrum for Dupuit aquifer Dirichlet
#------------------------------------------------------------------------------

# function of the analytical spectrum as a function of precipitation recharge
def func_dupuit_D(w, t_c):
    p = np.sqrt(1j*w*(t_c/(2.*np.pi)))
    
    return 1./(w*(t_c/(2.*np.pi)))*(abs(np.tanh(p))**2)

def func_dupuit_D_log(w, t_c): 
    return np.log10(func_dupuit_D(10**w, t_c))


def analytical_spectrum_dupuit_D(fr, w_QB, Pxx_QB):
    # spectrum baseflow
    (w_QB_log, Pxx_QB_log) = log_scale_spectrum(w_QB, Pxx_QB)
    
    popt, pcov = curve_fit(func_dupuit_D_log, w_QB_log, Pxx_QB_log, bounds=(0, np.inf), maxfev=1000)
    # Paramters fitted
    t_c = popt[0]

    return {'t_c':t_c}

#------------------------------------------------------------------------------
# Analytical spectrum for Dupuit aquifer Cauchy
#------------------------------------------------------------------------------

# function of the analytical spectrum as a function of precipitation recharge
def func_dupuit_C(w, t_c):
    p = np.sqrt(1j*w*(t_c/(2.*np.pi)))
    
    return 1./(w**2*(t_c/(2.*np.pi))**2)*(abs(np.tanh(p)/(np.tanh(p)+1/p))**2)

def func_dupuit_C_log(w, t_c): 
    return np.log10(func_dupuit_C(10**w, t_c))


def analytical_spectrum_dupuit_C(fr, w_QB, Pxx_QB):
    # spectrum baseflow
    (w_QB_log, Pxx_QB_log) = log_scale_spectrum(w_QB, Pxx_QB)
    
    popt, pcov = curve_fit(func_dupuit_C_log, w_QB_log, Pxx_QB_log, bounds=(0, np.inf), maxfev=1000)
    # Paramters fitted
    t_c = popt[0]

    return {'t_c':t_c}


#------------------------------------------------------------------------------
# Analytical spectrum for exponential spectrum for recharge
#------------------------------------------------------------------------------

# function of the analytical spectrum as a function of precipitation recharge
def func_exp(x, t_c, tau_R, beta):
    return np.log10((beta**2)/(1+tau_R**2*10**(2*x))/(1+t_c**2*10**(2*x)))

def analytical_spectrum_exp(fr, w_QB, Pxx_QB):
    # spectrum baseflow
    (w_QB_log, Pxx_QB_log) = log_scale_spectrum(w_QB, Pxx_QB)
    
    popt, pcov = curve_fit(func_exp, w_QB_log, Pxx_QB_log, p0=(1, 1, Pxx_QB_log[0]), bounds=(0, np.inf), maxfev=3000)

    t_c = max(popt[1],popt[0])*2.*np.pi    
    tau_R = min(popt[1],popt[0])*2.*np.pi
    beta = popt[2]

    # Calculate the goodness of fit    
    R = Rsquared(Pxx_QB_log,func_exp(w_QB_log, *popt))

    return {'analytical_spectrum': np.array([w_QB_log,func_exp(w_QB_log, *popt)]), 'Parameters':{'t_c':t_c, 'tau_R':tau_R, 'beta':beta, 'R':R}}

#------------------------------------------------------------------------------
# Analytical corrected spectrum
#------------------------------------------------------------------------------
    

def func_spectrum_corrected(x, t_c, beta, alpha):
    return np.log10((beta**2)/ (1  + t_c**2*10**(2*x))**alpha)

def analytical_spectrum_corrected(fr, w_QB, Pxx_QB):
    # spectrum baseflow
    (w_QB_log, Pxx_QB_log) = log_scale_spectrum(w_QB, Pxx_QB)
    
    popt, pcov = curve_fit(func_spectrum_corrected, w_QB_log, Pxx_QB_log, p0=(1000, np.sqrt(Pxx_QB_log[0]), 0), bounds=(0, np.inf), maxfev=1000)
    
    t_c = 2.*np.pi*popt[0]
    beta = popt[1]
    alpha = popt[2]

    # Calculate the goodness of fit    
    R = Rsquared(Pxx_QB_log,func_spectrum_corrected(w_QB_log, *popt))

    return {'analytical_spectrum': np.array([w_QB_log, func_spectrum_corrected(w_QB_log, *popt)]), 'Parameters':{'t_c':t_c, 'alpha':alpha, 'beta':beta, 'R':R}}