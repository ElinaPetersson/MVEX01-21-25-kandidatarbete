"""
Created on Thu May 13 18:48:04 2021.
@author: Jonas Lauri.
This file contains the ABC and MSE simulations for the Markov chain (CTMC) and the stochastic model (SDE).
"""

import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import random as random
import csv
import time

# %% The models. 
## The SIR model differential equations.
def SIR_data_Beta(y, t, N, beta_1, beta_2, gamma, f):
    S, I, R = y
    dSdt = -dataBeta(beta_1, beta_2, f, t) * S * I / N
    dIdt = dataBeta(beta_1, beta_2, f, t) * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

# Beta using mobility data.
def dataBeta(beta_1, beta_2, f, t):
    beta = beta_1 + beta_2 * f(t)
    return beta

## CTMC SIR model. 
def CTMC(N, beta, gamma, t_end, I0, R0):
    # Initialize the arrays.
    S_array = [0]
    I_array = [0]
    R_array = [0]
    t_array = [0]

    # Some starting values. 
    j = 0
    I_array[0], R_array[0], t_array[0], = I0, R0, 0
    S_array[0] = N - I_array[0] - R_array[0]

    # If I = 0 then the pandemic ends. 
    while I_array[j] > 0 and t_end > t_array[j]:
        a = beta * S_array[j] * I_array[j] / N
        b = gamma * I_array[j]
        # Compute the probability of a new infection.
        p_1 = a / (a + b)
        # 2 random numbers in the interval [0,1].
        u_1 = random.uniform(0,1)
        u_2 = random.uniform(0,1)
        if u_1 <= p_1:
            S_array.append(S_array[j] - 1)
            I_array.append(I_array[j] + 1)
            R_array.append(R_array[j])
        else:
            S_array.append(S_array[j])
            I_array.append(I_array[j] - 1)
            R_array.append(R_array[j] + 1)
        # Time of next event.
        t_array.append(t_array[j] - np.log(u_2) / (a + b))
        j = j + 1
    return S_array, I_array, R_array, t_array

## The SDE SIR-model, Euler-Maruyama method.
def SDE_EM(N, beta, gamma, dt, T, I0, R0):
    
    n = int(T / dt)  # Number of time steps.
    t = np.linspace(0, T, n)  # Vector of times.
    sqrtdt = np.sqrt(dt)

    # Initialize the arrays. 
    S = np.zeros(n)
    I = np.zeros(n)
    R = np.zeros(n)

    # Initial values.
    I[0], R[0] = I0, R0
    S[0] = N - I0 - R0

    # SDE definition
    # dS(t) = -[beta * S * I / N]dt - sqrt(beta * S * I / N)dW_1(t)
    # dI(t) = [beta * S * I / N - gamma * I]dt + sqrt(beta * S * I / N)dW_1(t) - sqrt(gamma * I)dW_2(t)
    # R = N - S - I

    # Euler-Maruyama method
    # dX(t) = f(X(t),t)dt + G(X(t),t)dW_1(t) ->
    # X(t + dt) = X(t) + f(X(t),t)dt + G(X(t),t)*eta*sqrt(dt)
    # where eta is a standard normal random number.

    # S(t + dt) = S(t) - [beta * S * I / N]dt - sqrt(beta * S * I / N)*eta*sqrt(dt)
    # I(t + dt) = I(t) + [beta * S * I / N - gamma * I]dt + sqrt(beta * S * I / N)eta*sqrt(dt) - sqrt(gamma * I)eta*sqrt(dt)
    # R(t + dt) = N - S(t + dt) - I(t + dt)

    for i in range(n - 1):
        eta_1, eta_2 = np.random.randn(), np.random.randn()
        S[i + 1] = S[i] - (beta * S[i] * I[i] / N)*dt + np.sqrt(beta * S[i] * I[i] / N)*eta_1*sqrtdt
        I[i + 1] = I[i] + (beta * S[i] * I[i] / N - gamma * I[i])*dt + np.sqrt(beta * S[i] * I[i] / N)*eta_1*sqrtdt - np.sqrt(gamma * I[i])*eta_2*sqrtdt
        R[i + 1] = N - S[i + 1] - I[i + 1]
    
    return S, I, R, t

# %% The data. 
# Hospitalization data. 
tsv_file = open("Covid19_sjukhus.tsv")
read_tsv = csv.reader(tsv_file, delimiter="\t")

# Extract the data to arrays.
week = []
currentHospital = []
cumulativeHospital = []
index = 0
# Extract the data to arrays.
for row in read_tsv:
    week.append(index)
    currentHospital.append(row[1])
    cumulativeHospital.append(row[2])
    index = index + 1
tsv_file.close()
# Want floats, not strings.
currentHospitalFloat = [float(currentHospital[i]) for i in range(len(currentHospital))]
cumulativeHospitalFloat = [float(cumulativeHospital[i]) for i in range(len(cumulativeHospital))]

# Turn weeks into days. 
days = np.linspace(0, 7*(len(week)-1), len(week))

# Interpolate the data. 
h_cur = interp1d(days, currentHospitalFloat, kind='linear', bounds_error = False, fill_value = 0)
h_cum = interp1d(days, cumulativeHospitalFloat, kind='linear', bounds_error = False, fill_value = 0) 

# Plot the current hospitalizations. 
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(days, h_cur(days))
ax.set_ylabel('Inläggningar')
ax.set_xlabel('Dagar')
plt.title('Nuvarande inläggningar', fontsize = 13)
plt.show()

# Plot the cumulative hospitalizations. 
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(days, h_cum(days))
ax.set_ylabel('Inläggningar')
ax.set_xlabel('Dagar')
plt.title('Kumulativa inläggningar', fontsize = 13)
plt.show()

## The mobility data.
tsv_file_2 = open("mobility_data.tsv")
read_tsv_2 = csv.reader(tsv_file_2, delimiter="\t")

# The same process for the other data set.
week = []
mobilityList = []
# Extract the data to arrays.
for row in read_tsv_2:
    week.append(row[0])
    mobilityList.append(row[2])
tsv_file_2.close()

# Fix the data.
mobilityData = [float(mobilityList[i]) for i in range(len(mobilityList))]
for i in range(len(mobilityList)):
    mobilityData[i] = 1 + mobilityData[i] / 100 
    i = i + 1
    
# Turn weeks into days. 
daysMob = np.linspace(0, 7*(len(week)-1), len(week))

# Interpolate the data. 
mobInterpol = interp1d(daysMob, mobilityData, kind='linear', bounds_error = False, fill_value = 0)

# Plot the mobility data. 
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(daysMob, mobInterpol(daysMob))
ax.set_ylabel('Mobilitetsindex')
ax.set_xlabel('Dagar')
plt.title('Mobilitetsdatan', fontsize = 13)
plt.show()

# %% ABC on the CTMC model.

# The ABC rejection algorithm:
# 1. Compute summary statistic my from observational data.
# 2. Given a certain model, perform n simulations, 
#    each with a parameter drawn from the prior the posterior distribution. 
# 3. Compute summary statistic my_i for each simulation. 
# 4. Based on a distance p(:,:) and a tolerance epsilon, 
#    decide for each simulation whether its summary statistic is suffiently close to that of the observed data.
# 5. Approximate the posterior distribution of theta from the distribution 
#    of parameter values theta_i associated with accepted simulations.

# The number of trials. 
n = 100000

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
scaling = 1/1000 # Simulation of N very large is not fun. 
N, gamma = 10000000 * scaling, 1./14
hospitalPar = 0.05 # How many infected people ends up needing hospital care?

# The maximum time duration of the outbreak.
t_end = 500

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0

# Save the values of the betas. 
betas = []

# Arbitrary limit on MSE. 
# MSE with only zeroes and scaling, 1/1000, is 1.4020050961665422.
mse_max = 1.25

# Prior bounds.
lowbnd, highbnd = 0.01, 1

# Returns MSE for given interpolated infected and interpolated data. 
# We have to use interpolated as the time is different each simulation. 
def MSE_interpol(days, data, I, scaling, hospitalPar):
    obj, i = 0, 0
    while i < len(days):
        obj = obj + (data(days[i])*scaling - I(days[i])*hospitalPar)**2
        i = i + 1
    return obj / i

# ABC algoritm.
now = time.time() # Timer.
for i in range(n):
    # Our prior. 
    theta = random.uniform(lowbnd, highbnd)
    
    # A trial run of CTMC.
    S_array, I_array, R_array, t_array = CTMC(N, theta, gamma, t_end, I0, R0)
            
    # Interpolate the I_array.
    f2 = interp1d(t_array, I_array, kind='linear', bounds_error = False, fill_value = 0)
        
    mse = MSE_interpol(days, h_cur, f2, scaling, hospitalPar)
        
    if mse < mse_max:
        betas.append(theta)
        # print(mse)

elapsed = time.time() - now
print(elapsed)

# Plot the resulting betas. 
plt.clf()
plt.hist(betas, bins=50)
plt.ylabel('Antal')
plt.xlabel('Beta')
plt.title('ABC för Markovkedjemodellen')
plt.plot([], [], ' ', label = "n = " + str(n))
plt.plot([], [], ' ', label = "mse_max = " + str(mse_max))
legend = plt.legend()
legend.get_frame().set_alpha(0.9)
plt.show()

# %% ABC and the SDE-model.

# The number of trials. 
n = 10000

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
scaling = 1
N, gamma = 10000000 * scaling, 1./14
hospitalPar = 0.05 # How many infected people ends up needing hospital care?

dt = 1  # Time step.
T = 385  # Total time.

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0

# Save the values of the betas. 
betas = []

# Arbitrary limit on MSE. 
# # MSE with only zeroes and scaling, 1/1000, is 1443883.1278705022.
mse_max = 1440000

lowbnd, highbnd = 0.01, 1

# A grid of time points.
t = np.linspace(0, 7*len(days) - 1)

# Returns MSE for given interpolated data.
def MSE(days, t, data, I, scaling, hospitalPar):
    obj, i = 0, 0
    while i < len(t):
        obj = obj + (data(t[i])*scaling - I[i]*hospitalPar)**2
        i = i + 1
    return obj / i

# ABC algoritm.
now = time.time() # Timer.
for i in range(n - 1):
    # Our prior. 
    theta = random.uniform(lowbnd, highbnd)
    
    # A trial run of SDE.
    S, I, R, t = SDE_EM(N, theta, gamma, dt, T, I0, R0)
    
    # Replace not a number with 0.
    I[np.isnan(I)] = 0
    
    mse = MSE(days, t, h_cur, I, scaling, hospitalPar)
        
    if mse < mse_max:
        betas.append(theta)
        print(mse)

elapsed = time.time() - now
print(elapsed)

# Plot the resulting betas. 
plt.clf()
plt.hist(betas, bins=50)
plt.ylabel('Antal')
plt.xlabel('Beta')
plt.title('ABC för SDE:n')
plt.plot([], [], ' ', label = "n = " + str(n))
plt.plot([], [], ' ', label = "mse_max = " + str(mse_max))
legend = plt.legend()
legend.get_frame().set_alpha(0.9)
plt.show()

# %% ABC and the SIR-model.

# The number of trials. 
n = 10000

# Total population, N and mean recovery rate, gamma, (in 1/days).
scaling = 1
N, gamma = 10000000, 1./14
hospitalPar = 0.05 # How many infected people ends up needing hospital care?

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
        
# Initial conditions vector.
y0 = S0, I0, R0

# Save the values of the betas. 
betas_1, betas_2 = [], []

# Arbitrary limit on MSE. 
mse_max = 1250000

lowbnd, highbnd = 0, 0.3

# A grid of time points.
t = np.linspace(0, 7*len(days) - 1)

now = time.time() # Timer.
for i in range(n):
    # Our priors. 
    theta_1 = random.uniform(lowbnd, highbnd)
    theta_2 = random.uniform(lowbnd, highbnd)
    
    ret = odeint(SIR_data_Beta, y0, t, args=(N, theta_1, theta_2, gamma, mobInterpol))
    S, I, R = ret.T
            
    mse = MSE(days, t, h_cur, I, scaling, hospitalPar)
        
    if mse < mse_max:
        betas_1.append(theta_1)
        betas_2.append(theta_2)
        # print(mse, theta_1, theta_2)

elapsed = time.time() - now
print(elapsed)

# Plot the resulting betas.
heatmap, xedges, yedges = np.histogram2d(betas_1, betas_2, bins = 20)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower')
plt.ylabel('beta_2')
plt.xlabel('beta_1')
plt.title('ABC för SIR-modellen med mobilitetsdata')
plt.plot([], [], ' ', label = "n = " + str(n))
plt.plot([], [], ' ', label = "mse_max = " + str(mse_max))
cbar = plt.colorbar()
cbar.ax.set_ylabel('')
legend = plt.legend()
legend.get_frame().set_alpha(0.5)
plt.show()

# %% Approximation of how many pandemics die early. 

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N_markov, N_sde, beta_markov, beta_sde, gamma = 1000, 1000, 0.1, 0.1, 1./14

# The time duration of the outbreak.
t_end = 500

dt = 1 # Time step.
T = 385 # Total time.

# Initial number of infected and recovered individuals, I0 and R0.
I0 = [1, 2, 3, 4, 5, 7, 10, 20]
R0 = 0

# Define what "dies early" means. 
largest_infected = [2, 4, 6, 8, 10, 14, 20, 40]

# Initialize the arrays.
propDeadCTMC = []
propDeadSDE = []

# Repeat the simulation n times. 
n = 1000

# Loop once for each I0.
now = time.time() # Timer.
for i in range(len(I0)):

    counter_ctmc, counter_sde = 0, 0
    
    for j in range(n - 1):
        S_array, I_array, R_array, t_array = CTMC(N_markov, beta_markov, gamma, t_end, I0[i], R0)
        S, I, R, t = SDE_EM(N_sde, beta_sde, gamma, dt, T, I0[i], R0)
        I_array_max = max(I_array)
        I_max = max(I)
        if largest_infected[i] > I_array_max:
            counter_ctmc = counter_ctmc + 1
        if largest_infected[i] > I_max:
            counter_sde = counter_sde + 1
            
    propDeadCTMC.append(counter_ctmc / n)
    propDeadSDE.append(counter_sde / n)

elapsed = time.time() - now
print(elapsed)

# Plot for proportion of dead pandemics with CTMC and SDE.
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(I0, propDeadCTMC, 'r', alpha=0.5, lw=2, label='CTMC')
ax.plot(I0, propDeadSDE, 'g', alpha=0.5, lw=2, label='SDE')
ax.set_xlabel('I0')
ax.set_ylabel('Sannolikhet')
plt.title('Sannolikheten att pandemin dör tidigt')
ax.plot([], [], ' ', label = "n = " + str(n))
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
plt.show()
