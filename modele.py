import numpy as np
from copy import copy
import matplotlib.pyplot as plt

class Strain:
    def __init__(self, a, n = "Strain basic"):
        # - Strain parameters --------------------------
        # Investement in aggregation
        self.alpha = a
        
        # Division rate --------------------------------
        self.c = 0.173

        # Scale factor of resources per division with cell size
        self.Q = 1

        # - Other attribute ----------------------------
        self.name = n
        self.spores = 0
        # Initial population
        self.vg = 10**8
        
    def __str__(self):
        return "Name : " + self.name + ", alpha = " + str(self.alpha) + ", Population X = " + str(self.vg + self.spores)
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

# - Global parameters - move this ---------------------

data = []  #stockage à chaque temps dt
R = R0  #ressources (va évoluer dans le temps)
strains = []
Time_passed = 0       
tau_counter = 0

# Normaly 2*(10**8)
total_time = 2*(10**3)

# -----------------------------------------------------

def initialisation(n_genotypes):
    """ Basic Simulation
    a = Strain(1)
    strains.append(a)
    """

    for i in range(n_genotypes):
        a = Strain(np.random.random(), "Strain " + str(i))
        strains.append(a)

def main_step(n, dt, n_starv):
    global Time_passed, R, tau_counter
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
    global Time_passed

    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    data.append(state)

    for s in strains: 
        #Calcul des conséquences de l'arrivée de la starvation
        s.spores = s.alpha * spore_stalk * s.vg
        s.vg = (1 - s.alpha) * s.vg

        #Evolution durant starvation
        s.spores -= delta * s.spores * t
        s.vg = survival_vg(t, s) * s.vg

def survival_vg(time_starvation, souche):
    global T_sur
    beta = 3.1 - (4 * souche.c)
    return (np.exp(-(mu*time_starvation)**(beta)) - np.exp(-(mu*T_sur)**(beta)))/(1 - np.exp(-(mu*T_sur)**(beta)))

def rich_step(dt):
    global Time_passed, R, tau_counter
    
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

def get_winner(d):
    max_pop = -1
    winner = None
    for s in d[-1]["strains"]:
        if s.vg + s.spores > max_pop:
            winner = s
            max_pop = s.vg + s.spores
    return winner


def plot_pop(d):
    strains_list = d[0]["strains"]
    t = []
    
    for i in range(len(d)):
        t.append(d[i]["Time_passed"])
    
    for k in range(len(strains_list)):
        arr = []
        for i in range(len(d)):
            arr.append(d[i]["strains"][k].vg + d[i]["strains"][k].spores)
        plt.plot(t, arr, label = strains_list[k].name)
    plt.show()

def plot_vgspore(d):
    strains_list = d[0]["strains"]
    t = []
    
    for i in range(len(d)):
        t.append(d[i]["Time_passed"])
    
    for k in range(len(strains_list)):
        arr = [[], []]
        for i in range(len(d)):
            arr[0].append(d[i]["strains"][k].vg)
            arr[1].append(d[i]["strains"][k].spores)
        plt.plot(t, arr[0], label = strains_list[k].name)
        plt.plot(t, arr[1], label = strains_list[k].name)
    plt.show()

def simulate_environements(max_starv_time, dt, n_genotypes = 1001, Tt = 2*(10**8), precision = 0.5):
    global data, R, strains, Time_passed, total_time, tau_counter
    # - Simulation
    total_time = Tt
    
    winners_alpha = []
    t = []
    
    for i in range(0, max_starv_time, dt):
        data = []
        R = R0
        strains = []
        Time_passed = 0       
        tau_counter = 0

        initialisation(n_genotypes)
        print("Simulation started ---")
        main_step(total_time, precision, i/precision)
        winners_alpha.append(get_winner(data).alpha)
        t.append(i)
        print("Currently on simulation: " + str(int(i/dt)) + "/" + str(int(max_starv_time/dt)))
    
    # - Graph
    plt.plot(t, winners_alpha, label = "alpha of winners")
    plt.show()

        
"""
initialisation() #Main
main_step(total_time, 1, 200)
print(get_winner(data))
"""
simulate_environements(400, 10)
print("helloworld")