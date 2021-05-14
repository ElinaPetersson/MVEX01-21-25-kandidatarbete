"""
Created on Sun Feb 28 13:26:47 2021
@author: Jonas Lauri, with the code for the basic SIR model being written by
https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/.
Additionally, the psuedocode for the CTMC SIR model is found here
https://mds.marshall.edu/cgi/viewcontent.cgi?article=2312&context=etd
Everything else is my own work. 
All parameters, files and models should be executable from this file. 
Add exposed as being able to infect?
Two exposed fack with one high beta and one with a low beta?
Vid smitta finns en viss sannolikhet att gå till varje fack?
Hur bestäms sannolikheten i så fall? Data på superspridare?
plt.savefig('destination_path.eps', format='eps')
Hur vill man göra med de fall då sjukdomen dör?
Fixa beta med varierande data. 
Finns tester för att undersöka posteriors fördelning, QQ etc. kolmogrov/smirnov
Scipy probablity. 
Måste se till MC anpassas till datan bättre. Tiden är ju olika varje gång!!!
MSE med bara 0: 2392.8736263736264!
Solve IVP i scipy.integrate, snabbare?
Göra sökandet mer systematiskt?
Behöver fixa så att beta beroende på data fungerar. 
Har ändrat gamma till 14 dagar. 
# Borde kunna fixa så att cumulativeHospitalFloat. 
Kanske? Bättre att använda currentHospitalFloat med hospitalRate. 
mse med bara 0: 1417291.2886276238.
"""
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
# from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import random as random
import csv
import time

# The basic SIR model. 
# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta, gamma = 1000, 0.2, 1./14

# A grid of time points (in days).
t = np.linspace(0, 160, 160)

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
        
# Initial conditions vector.
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
from Models import SIR
ret = odeint(SIR, y0, t, args=(N, beta, gamma))
# sol = solve_ivp(SIR, [0, 160], [990, 10, 0], args=(N, beta, gamma))
S, I, R = ret.T

from Plots import plot3
plot3(t, N, S, I, R)

# %% Markov Chain Epidemic Model.

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta, gamma = 1000, 0.2, 1./10

# The time duration of the outbreak.
t_end = 500

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0

from Models import CTMC
S_array, I_array, R_array, t_array = CTMC(N, beta, gamma, t_end, I0, R0)

from Plots import plot1
plot1(t_array, I_array)

# %% Compartmental model

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N_1, N_2, beta_1, beta_2, gamma_1, gamma_2 = 500, 300, 0.2, 0.1, 1./10, 1./15

# A grid of time points (in days).
t = np.linspace(0, 160, 160)

# Initial number of infected and recovered individuals, I0 and R0.
I1_0, I2_0, R_0 = 5, 5, 0
# Everyone else, S0, is susceptible to infection initially.
S1_0 = N_1 - I1_0
S2_0 = N_2 - I2_0

# The total population. 
N = N_1 + N_2

# Initial conditions vector.
y0 = S1_0, S2_0, I1_0, I2_0, R0

# Integrate the SIR equations over the time grid, t.
from Models import SIR2
ret = odeint(SIR2, y0, t, args=(N, beta_1, beta_2, gamma_1, gamma_2))
S1, S2, I1, I2, R = ret.T

from Plots import plot5
plot5(t, N, S1, S2, I1, I2, R)

# %% Time dependent beta. 

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta_0, gamma = 1000, 0.2, 1./10

# A grid of time points (in days).
t = np.linspace(0, 300, 300)

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
        
# Initial conditions vector.
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
from Models import SIR_time_Beta
ret = odeint(SIR_time_Beta, y0, t, args=(N, beta_0, gamma))
S, I, R = ret.T

plot3(t, N, S, I, R)


# %% The SEIR model, with constant incubation time. 

# Total population, N, contact rate, beta, mean recovery rate, gamma, (in 1/days) and incubation time, a. 
N, beta, gamma, a = 10000000, 0.2, 1./14, 7

# A grid of time points (in days).
t = np.linspace(0, 160, 160)

# Initial number of infected and recovered individuals, I0 and R0.
I0, E0, R0 = 2, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0

# Initial conditions vector.
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
from Models import SEIR
ret = odeint(SEIR, y0, t, args=(N, beta, gamma, a))
S, E, I, R = ret.T

# Plot the data for E(t).
plot1(t, E)

# %% The SDE SIR-model, Euler-Maruyama method.

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta, gamma = 1000, 0.2, 1./10

dt = 1  # Time step.
T = 300  # Total time.

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 2, 0

from Models import SDE_EM
S, I, R, t = SDE_EM(N, beta, gamma, dt, T, I0, R0)

# Plot the data on three separate curves for S(t), I(t) and R(t).
from Plots import plot3_SDE
plot3_SDE(t, N, S, I, R)

# %% Import the data from Folkhälsomyndigheten.
# Change below to current path to data!
tsv_file = open("AvlidnaPerDag_0309.tsv")
read_tsv = csv.reader(tsv_file, delimiter="\t")
# New data.
tsv_file_2 = open("Covid19_sjukhus.tsv")
read_tsv_2 = csv.reader(tsv_file_2, delimiter="\t")

# Extract the data to arrays.
date = []
deathsList = []
# Extract the data to arrays.
for row in read_tsv:
    date.append(row[0])
    deathsList.append(row[1])
tsv_file.close()
# Remove the deaths with unknown date. 
deathsList.pop(0)
deathsList = deathsList[:-1]
# Using map() to perform conversion.
deathsInt = list(map(int, deathsList)) 
# Normalize the dates. 
dates = np.linspace(0, len(date) - 3, len(date) - 2)

# Extract the data to arrays for the second data set.
week = []
currentHospital = []
cumulativeHospital = []
index = 0
# Extract the data to arrays.
for row in read_tsv_2:
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
h_cul = interp1d(days, cumulativeHospitalFloat, kind='linear', bounds_error = False, fill_value = 0) 

# Plot the result. 
from Plots import plot_data
plot_data(dates, deathsInt)
plot_data(days, h_cur(days))
plot_data(days, h_cul(days))

# %% Beta estimation for SIR using Mean Squared Error (MSE).
# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta, gamma = 1000, 0.01, 1./10

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0

# A grid of time points (in days).
t = np.linspace(0, len(dates)+1, len(dates))

# Initial conditions vector.
y0 = S0, I0, R0

# MSE_SIR(beta, gamma, N, dates, t, deathsInt, I, S0, I0, R0).
from MSE import MSE_SIR
result = minimize(MSE_SIR, [beta], args=(gamma, N, dates, t, deathsInt, y0))
#result = leastsq(MSE_SIR, [beta], args=(gamma, N, dates, t, deathsInt, y0), full_output = True)
print(result)
best_beta = result.x
#best_beta = result[0][0]

# Plot of solution using best_beta.
ret = odeint(SIR, y0, t, args=(N, best_beta, gamma))
S, I, R = ret.T

# Plot the result. 
from Plots import plot_data_model
plot_data_model(dates, deathsInt, t, I)

# %% Partition the data to get better results. 

# The size of the partition. 
size = 175

# Partition of deaths. 
from MSE import partition
deathsPart = list(partition(deathsInt, size))
deathsShort = list(map(int, deathsPart[0]))

# Partition of dates.
datesPart = list(partition(dates, size))
datesShort = list(map(int, datesPart[0]))

# A grid of time points (in days).
t = np.linspace(0, len(datesShort)+1, len(datesShort))

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta, gamma = 1000, 0.2, 1./10

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0

# Initial conditions vector.
y0 = S0, I0, R0

result = minimize(MSE_SIR, [beta], args=(gamma, N, dates, t, deathsInt, y0))
#print(result)
best_beta = result.x

# Plot of solution using best_beta.
ret = odeint(SIR, y0, t, args=(N, best_beta, gamma))
S, I, R = ret.T

# Plot the result. 
plot_data_model(datesShort, deathsShort, t, I)

# %% Beta estimation for compartamental model using MSE.

# Subpopulations N_1, N_2, contact rates, beta_1, beta_2, and mean recovery rates, gamma_1, gamma_2, (in 1/days).
N_1, N_2, beta_1, beta_2, gamma_1, gamma_2 = 500, 300, 0.2, 0.1, 1./10, 1./15

# A grid of time points (in days).
#t = np.linspace(0, len(dates)+1, len(dates))
t = np.linspace(0, len(datesShort)+1, len(datesShort))

# Initial number of infected and recovered individuals, I0 and R0.
I1_0, I2_0, R_0 = 5, 5, 0
# Everyone else, S0, is susceptible to infection initially.
S1_0 = N_1 - I1_0
S2_0 = N_2 - I2_0

# The total population.
N = N_1 + N_2

# Initial conditions vector.
y0 = S1_0, S2_0, I1_0, I2_0, R_0

# MSE_SIR2(beta_1, beta_2, gamma_1, gamma_2, N, dates, t, deathsInt, S1_0, S2_0, I1_0, I2_0, R_0).
# First minimize in the beta_1 direction. 
from MSE import MSE_SIR2_beta_1
#result = minimize(MSE_SIR2_beta_1, beta_1, args=(beta_2, gamma_1, gamma_2, N, dates, t, deathsInt, y0))
result = minimize(MSE_SIR2_beta_1, beta_1, args=(beta_2, gamma_1, gamma_2, N, datesShort, t, deathsShort, y0))
#print(result)
best_beta_1 = result.x

# Then minimize in the beta_2 direction.
from MSE import MSE_SIR2_beta_2
#result = minimize(MSE_SIR2_beta_2, beta_2, args=(beta_1, gamma_1, gamma_2, N, dates, t, deathsInt, y0))
result = minimize(MSE_SIR2_beta_2, beta_2, args=(beta_1, gamma_1, gamma_2, N, datesShort, t, deathsShort, y0))
#print(result)
best_beta_2 = result.x

# Plot of solution using best_beta.
ret = odeint(SIR2, y0, t, args=(N, best_beta_1, best_beta_2, gamma_1, gamma_2))
S1, S2, I1, I2, R = ret.T

# Plot the result. 
from Plots import plot_data_model_SIR2
plot_data_model_SIR2(datesShort, deathsShort, t, I1, I2)

# %% Beta estimation for SIR with time dependent beta. 

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta_0, gamma = 1000, 0.2, 1./10

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0

# A grid of time points (in days).
t = np.linspace(0, len(dates)+1, len(dates))

# Initial conditions vector.
y0 = S0, I0, R0

# MSE_SIR(beta, gamma, N, dates, t, deathsInt, I, S0, I0, R0).
from MSE import MSE_SIR_time
result = minimize(MSE_SIR_time, [beta_0], args=(gamma, N, dates, t, deathsInt, y0))
print(result)
best_beta = result.x

# Plot of solution using best_beta.
ret = odeint(SIR_time_Beta, y0, t, args=(N, best_beta, gamma))
S, I, R = ret.T

plot3(t, N, S, I, R)

# Plot the result. 
plot_data_model(dates, deathsInt, t, I)

# %% Approximation of how many pandemics die early. 

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N_markov, N_sde, beta_markov, beta_sde, gamma = 1000, 1000, 0.1, 0.1, 1./14

# The time duration of the outbreak.
t_end = 500

dt = 1  # Time step.
T = 385  # Total time.

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
    #print("Proportion dead CTMC: " + str(propDeadCTMC))
    propDeadSDE.append(counter_sde / n)
    #print("Proportion dead SDE: " + str(propDeadSDE))

from Plots import plot_data_branches
plot_data_branches(I0, propDeadCTMC, propDeadSDE, n)

elapsed = time.time() - now
print(elapsed)

# %% Beta estimation for SDE using MSE. Might be a lost cause. Commented.
"""
# This will done by averaging a number of simulations to counteract the randomness.
# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
N, beta, gamma = 1000, 0.2, 1./10
dt = 1  # Time step.
T = len(dates)  # Total time.
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 20, 0
# Repeat the simulation n times. 
n = 10
# A grid of time points (in days).
t = np.linspace(0, len(dates) - 1, len(dates))
from MSE import MSE_SDE
result = minimize(MSE_SDE, [beta], args=(gamma, N, dt, T, t, I0, R0, n, dates, deathsInt), options={'maxiter': 3, 'disp': True})
print(result)
best_beta = result.x
#S, I, R, t = SDE_EM(N, best_beta, gamma, dt, T, I0, R0)
#plot_data_model(dates, deathsInt, t, I)
"""

# %% SEAIR, 65+ compartment and vaccine usage. 

# Total population, N, contact rate, beta, mean recovery rate, gamma, (in 1/days) and incubation time, a.
N_1, N_2, beta_1, beta_2, gamma_1, gamma_2, a_1, a_2 = 500, 300, 0.2, 0.1, 1./10, 1./15, 2, 3

# Alpha is the rate at which asympomatic people spread the disease.
alpha_1, alpha_2 = 0.1, 0.1

# Sigma is the rate at which asymptomatic people move from asymtomatic to infected. 
sigma_1, sigma_2 = 0.1, 0.1

# Term for the speed of vaccination. 
q_1, q_2 = 0.01, 0.1

# A grid of time points (in days).
t = np.linspace(0, 160, 160)

# Initial number of exposed, infected, asymptomatic and recovered individuals, E0, I0, A0 and R0.
E1_0, E2_0, A1_0, A2_0, I1_0, I2_0, R_0 = 0, 0, 0, 0, 1, 1, 0
# Everyone else, S0, is susceptible to infection initially.
S1_0 = N_1 - I1_0 - E1_0
S2_0 = N_2 - I2_0 - E2_0

# The total population. 
N = N_1 + N_2

# Initial conditions vector.
y0 = S1_0, S2_0, E1_0, E2_0, A1_0, A2_0, I1_0, I2_0, R0

# Integrate the SIR equations over the time grid, t.
from Models import SEAIR2
ret = odeint(SEAIR2, y0, t, args=(N, beta_1, beta_2, gamma_1, gamma_2, a_1, a_2, alpha_1, alpha_2, sigma_1, sigma_2, q_1, q_2))
S1, S2, E1, E2, A1, A2, I1, I2, R = ret.T

from Plots import plotSEAIR1, plotSEAIR2
plotSEAIR1(t, N, S1, S2, R)
plotSEAIR2(t, N, E1, E2, A1, A2, I1, I2)

# %% Import the mobility data from Västtrafik.

# Change below to current path to data! 
tsv_file = open("mobility_data_old.tsv")
read_tsv = csv.reader(tsv_file, delimiter="\t")

# Change below to current path to data! 
tsv_file_2 = open("mobility_data.tsv")
read_csv = csv.reader(tsv_file_2, delimiter="\t")

# Extract the data to arrays.
week = []
mobilityList = []
# Extract the data to arrays.
for row in read_tsv:
    week.append(row[0])
    mobilityList.append(row[1])
tsv_file.close()

# Fix the data.
mobilityData = list(map(int, mobilityList)) 
for i in range(len(mobilityData)):
    mobilityData[i] = 1 + mobilityData[i] / 100
    i = i + 1
    
# Turn weeks into days. 
days = np.linspace(0, 7*(len(week)-1), len(week))

# Interpolate the data. 
f = interp1d(days, mobilityData, kind='linear', bounds_error = False, fill_value = 0)

# The same process for the other data set.
weekNew = []
mobilityListNew = []
# Extract the data to arrays.
for row in read_csv:
    weekNew.append(row[0])
    mobilityListNew.append(row[2])
tsv_file_2.close()

# Fix the data.
mobilityDataNew = [float(mobilityListNew[i]) for i in range(len(mobilityListNew))]
for i in range(len(mobilityListNew)):
    mobilityDataNew[i] = 1 + mobilityDataNew[i] / 100 # Är detta bästa sättet att använda datan?
    # mobilityDataNew[i] = 0.5 + mobilityDataNew[i] / 100
    i = i + 1
    
# Turn weeks into days. 
daysNew = np.linspace(0, 7*(len(weekNew)-1), len(weekNew))

# Interpolate the data. 
g = interp1d(daysNew, mobilityDataNew, kind='linear', bounds_error = False, fill_value = 0)

# Plot the result. 
from Plots import plot_data
plot_data(days, f(days))
plot_data(daysNew, g(daysNew))

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
n = 10000

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
scaling = 1/10000 # Simulation of N very large is not fun. 
N, gamma = 10000000 * scaling, 1./14
# N = 10000
hospitalPar = 0.05 # How many infected people ends up needing hospital care?

# The maximum time duration of the outbreak.
t_end = 500

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0

# Save the values of the betas. 
betas = []
# thetas = []

# Arbitrary limit on MSE. 
mse_max = 5

# A grid of time points (in days).
t = np.linspace(0, 7*len(days) - 1, len(days))

# Prior bounds.
lowbnd, highbnd = 0.01, 1

from Models import CTMC
from MSE import MSE_interpol_2

now = time.time() # Timer.
for i in range(n):
    # Our prior. 
    theta = random.uniform(lowbnd, highbnd)
    # thetas.append(theta)
    
    # A trial run of CTMC.
    S_array, I_array, R_array, t_array = CTMC(N, theta, gamma, t_end, I0, R0)
    
    I_array_max = max(I_array)
    
    # Want to quickly remove dying outbreaks. Make sure this makes sense. 
    if I_array_max > 50:
            
        # Interpolate the I_array.
        f2 = interp1d(t_array, I_array, kind='linear', bounds_error = False, fill_value = 0)
        
        mse = MSE_interpol_2(dates, t, h_cur, f2, scaling, hospitalPar)
        
        if mse < mse_max:
            betas.append(theta)
            #print(mse)

elapsed = time.time() - now
print(elapsed)

# Plot the resulting betas. 
plt.clf()
# plt.hist(thetas, bins=50)
plt.hist(betas, bins=50)
plt.ylabel('Count')
plt.xlabel('Beta')
plt.title('ABC for the CTMC model')
plt.plot([], [], ' ', label = "n = " + str(n))
plt.plot([], [], ' ', label = "mse_max = " + str(mse_max))
legend = plt.legend()
legend.get_frame().set_alpha(0.5)
plt.show()

# from Plots import plot1
# plot_data(days, scaling*h_cur(days))
# plot1(t_array, I_array)

# %% ABC and the SDE-model.

# The number of trials. 
n = 100000

# Total population, N, contact rate, beta, and mean recovery rate, gamma, (in 1/days).
scaling = 1 # Simulation of N very large is not fun. 
N, gamma = 10000000 * scaling, 1./14
hospitalPar = 0.05 # How many infected people ends up needing hospital care?

dt = 1  # Time step.
T = 385  # Total time.

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0

# Save the values of the betas. 
betas = []

# Arbitrary limit on MSE. 
mse_max = 1500000 * scaling

lowbnd, highbnd = 0.01, 1

# A grid of time points (in days).
t = np.linspace(0, 7*len(days) - 1, len(days))

from Models import SDE_EM
from MSE import MSE_2

now = time.time() # Timer.
for i in range(n - 1):
    # Our prior. 
    theta = random.uniform(lowbnd, highbnd)
    
    # A trial run of SDE.
    S, I, R, t = SDE_EM(N, theta, gamma, dt, T, I0, R0)
    
    # Replace not a number with 0.
    I[np.isnan(I)] = 0
    
    I_max = max(I)
    
    # Want to remove quickly dying outbreaks. 
    if I_max > 50:
            
        mse = MSE_2(dates, t, h_cur, I, scaling, hospitalPar)
        
        if mse < mse_max:
            betas.append(theta)
            # print(mse)

elapsed = time.time() - now
print(elapsed)

# Plot the resulting betas. 
plt.clf()
plt.hist(betas, bins=50)
plt.ylabel('Count')
plt.xlabel('Betas')
plt.title('ABC for the SDE model')
plt.plot([], [], ' ', label = "n = " + str(n))
plt.plot([], [], ' ', label = "mse_max = " + str(mse_max))
legend = plt.legend()
legend.get_frame().set_alpha(0.5)
plt.show()

from Plots import plot3_SDE
# plot3_SDE(t, N, S, I, R)

# %% ABC and the SIR-model.

# The number of trials. 
n = 10

# Total population, N and mean recovery rate, gamma, (in 1/days).
scaling = 1 # Simulation of N very large is not fun. 
# N, gamma = 10000000, 1./14
N, gamma = 10000, 1./14
hospitalPar = 0.05 # How many infected people ends up needing hospital care?

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 10, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
        
# Initial conditions vector.
y0 = S0, I0, R0

# Save the values of the betas. 
betas_1, betas_2 = [], []
S_saved = [] 

# Arbitrary limit on MSE. 
mse_max = 1250000
# mse_max = 1

lowbnd, highbnd = 0, 0.3

# A grid of time points (in days).
# t = np.linspace(0, 7*len(daysNew), len(daysNew)+1)
t = np.linspace(0, 7*25, 25+1)

from Models import SIR_data_Beta
from MSE import MSE_2

now = time.time() # Timer.
for i in range(n):
    # Our priors. 
    theta_1 = random.uniform(lowbnd, highbnd)
    theta_2 = random.uniform(lowbnd, highbnd)
    
    ret = odeint(SIR_data_Beta, y0, t, args=(N, theta_1, theta_2, gamma, g))
    S, I, R = ret.T
            
    mse = MSE_2(days, t, h_cur, I, scaling, hospitalPar)
    
    S_saved.append(S)
        
    if mse < mse_max:
        betas_1.append(theta_1)
        betas_2.append(theta_2)
        print(mse, theta_1, theta_2)

elapsed = time.time() - now
print(elapsed)

# Plot the resulting betas.
heatmap, xedges, yedges = np.histogram2d(betas_1, betas_2, bins = 20)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower')
plt.ylabel('beta_2')
plt.xlabel('beta_1')
# Create empty plot with blank marker containing the extra label.
plt.plot([], [], ' ', label = "n = " + str(n))
plt.plot([], [], ' ', label = "mse_max = " + str(mse_max))
legend = plt.legend()
legend.get_frame().set_alpha(0.5)
# plt.axis('square')
plt.show()

plt.hist2d(betas_1, betas_2, bins=20)
plt.ylabel('betas_2')
plt.xlabel('betas_1')
plt.title('ABC for the SIR model with mobility data')
plt.plot([], [], ' ', label = "n = " + str(n))
plt.plot([], [], ' ', label = "mse_max = " + str(mse_max))
cbar = plt.colorbar()
cbar.ax.set_ylabel('')
# plt.savefig('destination_path.eps', format='eps')
plt.show()

# Plot of the mean of beta_1 and beta_2. TODO!
# Används i en plot av Joakim.
beta_1, beta_2 = 0.1, 0.09

# Solve the SIR.
ret = odeint(SIR_data_Beta, y0, t, args=(N, beta_1, beta_2, gamma, g))
S, I, R = ret.T

from Plots import plot3
plot3(t, N, S, I, R)
