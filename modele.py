import numpy as np
from copy import copy
import matplotlib.pyplot as plt

class Strain:
    def __init__(self, a):
        # - Strain parameters --------------------------
        # Investement in aggregation
        self.alpha = a
        
        # Division rate --------------------------------
        self.c = 0.173

        # Scale factor of resources per division with cell size
        self.Q = 1

        # - Other attribute ----------------------------
        self.name = "Strain basic"
        self.spores = 0
        # Initial population
        self.vg = 10**8
        
    def __str__(self):
        return "Name : " + self.name + ", Population X = " + str(self.vg + self.spores)
    def __repr__(self):
        return self.__str__()

# - Copy strains without reference
def copyStrains(strains):
    resultat = []
    for s in strains:
        resultat.append(copy(s))
    return resultat

# - Simulation parameters (described in paper 108) -----
# Rate of decrease of the survival probability
mu = 2*(10**-3)

# Maximum lifetime of a non-aggregated starving cell --
T_sur = 200

# Spore germination time ------------------------------
tau = 4

# Spore mortality rate --------------------------------
delta = 2 * (10**(-4))

# Fraction of aggregators that become spores ----------
spore_stalk = 0.8

# Food Pusle Size -------------------------------------
R0 = 10**8

# Half-saturation constant of resources consumption ---
Rhalf = 0.1*R0

# Resources exhaustion threshold ----------------------
Rstar = 1

# - Global parameters ---------------------------------

data = []  #stockage à chaque temps dt
R = R0  #ressources (va évoluer dans le temps)
strains = []
Time_passed = 0       
tau_counter = 0
Time_rich = 0
Time_starvation = 0

total_time = 2*(10**3)

# -----------------------------------------------------

def initialisation():
    a = Strain(0)
    strains.append(a)

def main_step(n, dt, n_starv):
    global Time_passed, R, tau_counter, Time_rich, Time_starvation
    i = 0
    while i < n * dt:
        if R > Rstar:
            rich_step(dt)
            i += dt
        else:
            starvation_step(n_starv * dt)
            i += n_starv * dt
            R = 10**8
            tau_counter = 0

def starvation_step(t):
    global Time_passed, Time_starvation, Time_rich

    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    data.append(state)

    Time_rich = 0

    for s in strains: 
        #Calcul des conséquences de l'arrivée de la starvation
        s.spores = s.alpha * spore_stalk * s.vg
        s.vg = (1 - s.alpha) * s.vg

        #Evolution durant starvation
        s.spores -= delta * s.spores * t
        s.vg = survival_vg(t, s) * s.vg

    Time_starvation += t
    Time_passed += t

def survival_vg(time_starvation, souche):
    global T_sur
    beta = 3.1 - (4 * souche.c)
    return (np.exp(-(mu*time_starvation)**(beta)) - np.exp(-(mu*T_sur)**(beta)))/(1 - np.exp(-(mu*T_sur)**(beta)))

def rich_step(dt):
    global Time_passed, R, tau_counter, Time_rich, Time_starvation
    #print(Time_passed)

    Time_starvation = 0 #! ici ça doit bien être temps de starvation
    
    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    data.append(state)

    #spores germination
    if tau_counter >= tau and state['Time_passed']!=0:
        for s in strains:
            #On calcule le flux de nb de spores germinant en vg
            nu = 1.1 - 2*s.c
            s.vg += s.spores * (1 - nu)
            s.spores = 0


            ###Sens inverse :
            # flow = s.vg * s.alpha #!!! alpha bien par rapport à vg slmt ?
            # s.vg -= flow #On enlève la part correspondant aux cellules du sporocarpe, qui ne deviennent pas des endospores et ne se reproduisent pas (donc négligeables)
            # s.spores += flow * spocre_stalk

        tau_counter = -1
    
    elif tau_counter >= tau : #(and Time_passed == 0)
        print("pas Germination car situation initiale")

    #calcul
    kRGab = state["R"] / (state["R"] + Rhalf)
    sommeTit = 0

    #calcul - population
    for s in strains:
        X = s.vg + s.spores
        sommeTit += s.Q * s.c * X
        s.vg += s.c * kRGab * X * dt
        if s.spores > 0: s.spores -= delta * s.spores * dt

    # calcul - ressource
    R += - kRGab * sommeTit
    if R <= Rstar:
        R = Rstar

    Time_passed += dt
    if tau_counter >= 0: tau_counter += dt
    
    Time_rich += dt



initialisation() #Main
main_step(total_time, 1, 100)
Pop1R_plot = [[],[],[],[]]

def plot_popstrain(strain_name):
    strains_list = data[0]["strains"]
    index_strain = 0
    for j in range(len(strains_list)):
        if strains_list[j]==strain_name:
            index_strain = j
    t = []
    for i in range(len(data)):
        Pop1R_plot[0].append(data[i]["R"])
        Pop1R_plot[1].append(data[i]["strains"][index_strain].vg + data[i]["strains"][index_strain].spores)
        Pop1R_plot[2].append(data[i]["strains"][index_strain].vg)
        Pop1R_plot[3].append(data[i]["strains"][index_strain].spores)
        t.append(data[i]["Time_passed"])

    
    
    #plt.scatter(t, Pop1R_plot[0], label='Ressources')
    #plt.scatter(t, Pop1R_plot[1], label="Population")

    #plt.plot(t, Pop1R_plot[0], label="Ressources")
    #plt.plot(t, Pop1R_plot[1], label="Population")
    plt.plot(t, Pop1R_plot[2], label="Vege")
    plt.plot(t, Pop1R_plot[3], label="Spores")

    plt.show()
    plt.legend(['First line', 'Second line'])

plot_popstrain('Strain basic')

print("helloworld")