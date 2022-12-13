import numpy as np
import math
import json
from scipy.special  import erfinv
from copy import copy
import matplotlib.pyplot as plt

class Strain:
    def __init__(self, a, n = "Strain basic", r = 1):
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
        self.initial_ratio = r
    def normalize_pop(self, normalize):
        self.initial_ratio = self.initial_ratio/normalize
        self.vg = 10**8 * self.initial_ratio
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

data = []  # stockage à chaque temps dt
R = R0  # -- ressources (va évoluer dans le temps)
strains = []
Time_passed = 0       
tau_counter = 0

# Normaly 2*(10**8)
total_time = 2*(10**4)

# -----------------------------------------------------

def initialisation(n_genotypes):
    """
    Basic Simulation
    a = Strain(1)
    strains.append(a)
    """
    lognormal = lambda m, v, p : np.exp(m + np.sqrt(2 * v) * erfinv(2 * p - 1))
    logsum = 0
    
    for i in range(n_genotypes):
        x = lognormal(0, 1, np.random.rand())
        a = Strain(np.random.random(), "Strain " + str(i), x)
        strains.append(a)
        logsum += x
    
    for s in strains:
        s.normalize_pop(logsum)

def main_step(n, dt, n_starv, quant_func = lambda x, p : x):
    global Time_passed, R, tau_counter
    i = 0
    while i < n * dt:
        if R > Rstar:
            rich_step(dt)
            i += dt
        else:
            starv_time = quant_func(n_starv, np.random.rand())
            if starv_time > 0:
                starvation_step(starv_time)
            i += n_starv
            R = 10**8
            tau_counter = 0

def wetseveredry_step(n, dt, T_dry, lambda_wet, lambda_wet_function):
    global Time_passed, R, tau_counter
    T_wet = 8760 - T_dry
    dry_time = False
    
    i = 0
    
    while i < n * dt:
        if not dry_time:
            main_step(T_wet/dt, dt, lambda_wet, lambda_wet_function)
            dry_time = True
            
            i += T_wet
        
        else:
            starvation_step(T_dry)
            dry_time = False
            
            if T_wet > 1: R = 10**8
            tau_counter = 0

            i += T_dry
            
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
        
        if s.vg < 0: s.vg = 0
        if s.spores < 0: s.spores = 0
    Time_passed += t

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
        else: s.spores = 0

    # calcul - ressource
    R += -kRGab * sommeTit * dt
    if R <= Rstar:
        for i in range(len(strains)):
            strains[i].vg = state["strains"][i].vg
            strains[i].spores = state["strains"][i].spores
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

def get_winner_alpha_mean(d):
    avg = 0
    tpop = 0
    for s in d[-1]["strains"]:
        avg += s.alpha * (s.vg + s.spores)
        tpop += s.vg + s.spores
    return avg / tpop

def plot_pop(d, title="Evolution des population"):
    strains_list = d[0]["strains"]
    t = []
    
    for i in range(len(d)):
        t.append(d[i]["Time_passed"])
    
    for k in range(len(strains_list)):
        arr = []
        for i in range(len(d)):
            arr.append(d[i]["strains"][k].vg + d[i]["strains"][k].spores)
        plt.plot(t, arr, label = "α = {x:.2f}".format(x = strains_list[k].alpha))
    
    plt.title(title, fontsize=8)
    plt.xlabel("Temps (heure)")
    plt.ylabel("Taille de population, spores + végétatives")
    plt.legend(loc = "lower right", fontsize=6)
    plt.show()

def plot_vgspore(d, title = "None"):
    strains_list = d[0]["strains"]
    t = []
    
    for i in range(len(d)):
        t.append(d[i]["Time_passed"])
    
    for k in range(len(strains_list)):
        arr = [[], []]
        for i in range(len(d)):
            arr[0].append(d[i]["strains"][k].vg)
            arr[1].append(d[i]["strains"][k].spores)
        plt.plot(t, arr[0], label = "α = {x:.2f}, végétatives".format(x = strains_list[k].alpha))
        plt.plot(t, arr[1], label = "α = {x:.2f}, spores".format(x = strains_list[k].alpha))
    
    plt.title(title, fontsize=8)
    plt.xlabel("Temps (heure)")
    plt.ylabel("Taille de population")
    plt.legend(loc = "lower right")
    plt.show()

def simulate_environements(max_starv_time, dt, n_genotypes = 1001, Tt = 2*(10**2), precision = 0.1, function = lambda x, p : x, output_name = "output.txt"):
    global data, R, strains, Time_passed, total_time, tau_counter
    # - Simulation
    total_time = Tt
    
    winners_alpha = []
    t = []
    
    for i in range(1, max_starv_time, dt):
        data = []
        R = R0
        strains = []
        Time_passed = 0       
        tau_counter = 0

        initialisation(n_genotypes)
        print("Simulation started ---")
        main_step(total_time, precision, i, function)
        winners_alpha.append(get_winner_alpha_mean(data))
        t.append(i)
        print("Currently on simulation: " + str(int(i/dt) + 1) + "/" + str(int(max_starv_time/dt)))
    
    # - Graph
    plt.plot(t, winners_alpha, label = "alpha of winners")

    # - Save data just in case
    output = [t, winners_alpha]
    save_file_with_data(output_name, output)

def simulate_severe_dry_environements(max_wet_time, dt, T_dry, lambda_wet_function, n_genotypes = 1001, Tt = 2*(10**2), precision = 0.1, output_name = "output.txt"):
    global data, R, strains, Time_passed, total_time, tau_counter
    # - Simulation
    total_time = Tt
    
    winners_alpha = []
    t = []
    
    for i in range(1, max_wet_time, dt):
        data = []
        R = R0
        strains = []
        Time_passed = 0       
        tau_counter = 0

        initialisation(n_genotypes)
        print("Simulation started ---")
        #main_step(total_time, precision, i, function)
        wetseveredry_step(total_time, precision, T_dry, i, lambda_wet_function)
        winners_alpha.append(get_winner_alpha_mean(data))
        t.append(i)
        print("Currently on simulation: " + str(int(i/dt) + 1) + "/" + str(int(max_wet_time/dt)))
    
    # - Graph
    plt.plot(t, winners_alpha, label = "alpha of winners")

    # - Save data just in case
    output = [t, winners_alpha]
    save_file_with_data(output_name, output)

def save_file_with_data(n, d):
    f = open("./outputs/" + n, "w")
    f.write(json.dumps(d))
    f.close

determinist = lambda x, p : x
exponential = lambda x, p : -np.log(1 - p) * x
normal_20 = lambda x, p : x + 20 * np.sqrt(2) * erfinv(2 * p - 1)
normal_40 = lambda x, p : x + 40 * np.sqrt(2) * erfinv(2 * p - 1)
normal_60 = lambda x, p : x + 60 * np.sqrt(2) * erfinv(2 * p - 1)

data = []
R = R0
strains = []
Time_passed = 0       
tau_counter = 0

initialisation(10)
print("Simulation started ---")
#main_step(total_time, precision, i, function)
main_step(1*(10**5), 0.1, 200, determinist)
print(get_winner(data))
plot_pop(data, "Simulation, $T_{sur} =$ 200 heures (déterministe)")

data = []
R = R0
strains = []
Time_passed = 0       
tau_counter = 0

initialisation(10)
print("Simulation started ---")
#main_step(total_time, precision, i, function)
main_step(1*(10**5), 0.1, 200, exponential)
print(get_winner(data))
plot_pop(data, "Simulation, $T_{sur}$ suis une loi exponentielle de moyenne 200 heures.")

data = []
R = R0
strains = []
Time_passed = 0       
tau_counter = 0

initialisation(10)
print("Simulation started ---")
#main_step(total_time, precision, i, function)
main_step(1*(10**5), 0.1, 200, normal_20)
print(get_winner(data))
plot_pop(data, "Simulation, $T_{sur}$ suis une loi normale de moyenne 200 heures et d'écart-type 20.")

data = []
R = R0
strains = []
Time_passed = 0       
tau_counter = 0

initialisation(10)
print("Simulation started ---")
#main_step(total_time, precision, i, function)
main_step(1*(10**5), 0.1, 200, normal_60)
print(get_winner(data))
plot_pop(data, "Simulation, $T_{sur}$ suis une loi normale de moyenne 200 heures et d'écart-type 60.")
print("hello")
