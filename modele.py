import numpy as np
from copy import copy

class Strain:
    def __init__(self):
        self.name = "Strain basic"
        self.c = 0.173        #division rate
        self.Q = 1         #a voir, cmax/c
        self.X = 1000        #population
        self.alpha = 0.5
    def __str__(self):
        return "Name : " + self.name + ", Population X = " + str(self.X)
    def __repr__(self):
        return self.__str__()

# Copy strains without reference
def copyStrains(strains):
    resultat = []
    for s in strains:
        resultat.append(copy(s))
    return resultat

#global parameters
data = []
R0 = 10    #ressources
Rhalf = 0.1*R0
R = R0
Rstar = 1
strains = []
Time_passed = 0       
tau = 4       #hours
tau_counter = 0
delta = 2 * 10**(-4)    #spore mortality rate

def initialisation() :
    a = Strain()
    strains.append(a)

def rich_step(dt):
    global Time_passed, R, tau_counter
    #print(Time_passed)
    
    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    data.append(state)

    #spores germination
    if tau_counter >= tau:
        print("Germination")
        for s in strains:
            nu = 1.1 - 2*s.c
            s.X -= s.X * s.alpha * (1 - nu)
        tau_counter = -1

    #calcul
    kRGab = state["R"] / (state["R"] + Rhalf)
    sommeTit = 0

    #calcul - population
    for s in strains:
        sommeTit += s.Q * s.c * s.X
        s.X += s.c * kRGab * s.X * dt
        if tau_counter < tau: s.X -= delta * s.X * s.alpha * dt

    # calcul - ressource
    R += - kRGab * sommeTit
    if R <= Rstar:
        R = Rstar

    Time_passed += dt
    if tau_counter >= 0: tau_counter += dt

initialisation()
for i in range(1,600) :
    rich_step(1)

print(data)

print("helloworld")
