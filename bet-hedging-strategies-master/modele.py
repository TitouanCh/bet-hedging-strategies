import numpy as np
from copy import copy
import matplotlib.pyplot as plt

class Strain:
    def __init__(self):
        self.name = "Strain basic"
        self.c = 0.173        #division rate
        self.Q = 1         #a voir, cmax/c
        self.nbSpores = 0
        self.nbVG = 1000
        self.X = self.nbSpores + self.nbVG        #population
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
R0 = 10**8    #ressources
Rhalf = 0.1*R0
R = R0
Rstar = 1
strains = []
Time_passed = 0       
tau = 4       #hours
tau_counter = 0
delta = 2 * 10**(-4)    #spore mortality rate
spore_stalk = 0.8

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
    if tau_counter >= tau and state['Time_passed']!=0:
        print("Germination")
        for s in strains:
            nu = 1.1 - 2*s.c
            flow = s.nbVG * s.alpha * (1 - nu) #!!! alpha bien par rapport à VG slmt ?
            s.X -= s.X * s.alpha * (1 - nu)

        tau_counter = -1
    
    elif tau_counter >= tau :
        print("pas Germination car situation initiale")

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

initialisation() #Main
for i in range(1,600) :
    rich_step(1)
Pop1R_plot = [[],[]]

def plot_popstrain(strain_name):
    strains_list = data[0]["strains"]
    index_strain = 0
    for j in range(len(strains_list)):
        if strains_list[j]==strain_name:
            index_strain = j
    for i in range(len(data)):
        Pop1R_plot[0].append(data[i]["R"])
        Pop1R_plot[1].append(data[i]["strains"][index_strain].X)
    
    t = np.arange(len(Pop1R_plot[0]))
    
    #plt.scatter(t, Pop1R_plot[0], label='Ressources')
    #plt.scatter(t, Pop1R_plot[1], label="Population")

    plt.plot(t, Pop1R_plot[0])
    plt.plot(t, Pop1R_plot[1])

print(plot_popstrain('Strain basic'))

print(data)

print("helloworld")
