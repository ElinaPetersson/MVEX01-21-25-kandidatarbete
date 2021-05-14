'''
ABC-simulering för nätverksmodellen.
Input är modellen och output är en posteriorifördelning för parametrerna pSE och weightPar.

Filen reloadData måste köras först.
'''

from networkModel  import * 
import seaborn as sns
import numpy as np
import time
import matplotlib.pyplot as plt

# Beräknar avstånd mellan modellen och observerade värden.
def distance(y_model, y_observed):
    dis = 0
    for i in range(len(y_observed)):
        dis = dis + (y_model[i] - y_observed[i]) ** 2
    return dis

# Beräknar en avstånd mellan observerade värden och inga sjukhusinläggningar. Detta sätts till tolerans för ABC:n.
def tolCalculator():
    y_observed = model(0,0,False)['y_hospital_observed']
    y_model = [0 for i in range(len(y_observed))]
    tol = distance(y_model, y_observed)
    return tol

# ABC som sparar accepterade värden för parameterkombinationer.
def ABC(topPSE, topWeightPar, tol, simulations):
    acceptedPSE = []
    acceptedWeight = []
    for j in range(simulations):
        testPSE = random.uniform(0, topPSE)
        testWeightPar = random.uniform(0, topWeightPar)
        mod = model(testPSE, testWeightPar, False)
        if distance(mod['y_hospital_model'], mod['y_hospital_observed']) < tol:
            acceptedPSE.append(testPSE)
            acceptedWeight.append(testWeightPar)
        print(j)
    return [acceptedPSE,acceptedWeight]

# plottar upp det ABC resulterat i.        
def plotHeatMap(acceptedPSE,acceptedWeight):
    plt.hist2d(acceptedPSE,acceptedWeight,bins=20)
    plt.ylabel('weightPar')
    plt.xlabel('pSE')
    plt.title('ABC för nätverksmodellen', fontsize=13)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('')
    plt.show()
            
t0 = time.time()
tolFactor = 0.3 # nu är den 3 gånger bättre än inga sjukhusinläggningar alls
topPSE = 0.2 # övre gräns i apriorifördelningen
topWeightPar = 1 # övre gräns i apriorifördelningen
simulations = 10000

[acceptedPSE, acceptedWeight] = ABC(minPSE, minWeightPar, tolCalculator() * tolFactor, simulations)
plotHeatMap(acceptedPSE, acceptedWeight)


print('Tot tid: ', time.time() - t0)
