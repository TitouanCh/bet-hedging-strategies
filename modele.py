import numpy as np
from copy import copy
import matplotlib.pyplot as plt
import random

class Strain:
    def __init__(self):
        self.name = "Strain basic"
        self.c = 0.173        #division rate
        self.Q = 1         #a voir, cmax/c
        self.spores = 0
        self.vg = 10**8       #population
        self.alpha = 0.5
    def __str__(self):
        return "Name : " + self.name + ", Population X = " + str(self.vg + self.spores)
    def __repr__(self):
        return self.__str__()

# Copy strains without reference
def copyStrains(strains):
    resultat = []
    for s in strains:
        resultat.append(copy(s))
    return resultat

#global parameters
data = []  #stockage à chaque temps dt
R0 = 10**8    #ressources initiales
Rhalf = 0.1*R0  #moitié de saturation des ressources
R = R0  #ressources (va évoluer dans le temps)
Rstar = 1
strains = []
# new_strains = []
winners = [] #Les alphas des meilleures strains pour chaque sélection (si tout va bien, dans le même ordre que celui donné par les curseurs des autres paramètres)
Time_passed = 0       
tau = 4       #hours = spore germination time
tau_counter = 0
delta = 2 * 10**(-4)    #spore mortality rate
spore_stalk = 0.8
Time_rich = 0
Time_starvation = 0
T_sur = 200


def survival_vg(time_starvation, souche):
    global T_sur
    mu = 1.1 - (2 * souche.c)
    beta = 3.1 - (4 * souche.c)
    return (np.exp(-(mu*time_starvation)**(beta)) - np.exp(-(mu*T_sur)**(beta)))/(1 - np.exp(-(mu*T_sur)**(beta)))

""" def initialisation():
    a = Strain()
    strains.append(a) """

def Normale_density(x, mu, sigma):
    return (np.pi*sigma) * np.exp(-0.5*((x-mu)/sigma)**2)

def Normale_starvationtime_randomselection(parameters_mu_sigma, step=0.05):                  #(OLD : domain = borne sup (borne inf = 0), size ->~= precision)
    mu = parameters_mu_sigma[0]
    sigma = parameters_mu_sigma[1]
    
    x=0
    prob_table = [Normale_density(x, mu, sigma)]
    x+=1

    while prob_table[-1]>=0.001:
        prob_table.append(Normale_density(x, mu, sigma))
        x+=step

    i = random.randint(0, len(prob_table))
    return prob_table[i]

""" print(Normale_starvationtime_randomselection([0,1]))
tab = Normale_starvationtime_randomselection([0,1])
plt.plot(np.arange(len(tab)), tab)
plt.show() """

# def main_step(n, dt, n_starv): #
#     global Time_passed, R, tau_counter, Time_rich, Time_starvation
#     i = 0
#     while i < n * dt:
#         if R > Rstar:
#             rich_step(dt)
#             i += dt
#         else:
#             starvation_step(n_starv * dt)
#             i += n_starv * dt
#             R = 10**8
#             tau_counter = 0

def main_step(n, dt, parameters, density=Normale_starvationtime_randomselection, n_starv=None): #
    global Time_passed, R, tau_counter, Time_rich, Time_starvation
    i = 0
    while i < n * dt:
        if R > Rstar:
            rich_step(dt)
            i += dt
        else:
            #normale_starvation()
            n_starv = density(parameters)
            starvation_step(n_starv * dt) #Remplacer par fonction répartition - "the growth phase starts with the arrival of a pulse of ressources"
            i += n_starv * dt
            R = 10**8 #Pulse of resources
            tau_counter = 0

#def normale_starvation(sigma):
    #return un temps de starvation -> voir temps moyen de starvation

#    pass

def starvation_step(t):
    global Time_passed, Time_starvation, Time_rich

    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    ### data.append(state)

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



""" initialisation() #Main
main_step(400, 0.5, 4)
Pop1R_plot = [[],[],[],[]] """
Pop1R_plot = [[],[],[],[]]


def Exp_starvationtime(lambdaST, size, domain):
    pass


def Simulation_and_Selection(parameters, density=Normale_starvationtime_randomselection):   #Lance toute la simulation pour une densité de proba dépendant des params donnés,
                                                                                                     #et sélectionne le meilleur génotype
    global strains
    global winners
    
    graph_return = []
    #tmax = 2*108
    tmax = 10 ####
    dt = 1 #h

    main_step(tmax, dt, parameters, density)

    counting_survivors = []
    for s in range(len(strains)):
        counting_survivors.append(data[-1]["strains"][s].spores + data[-1]["strains"][s].vg)

    best = 0
    for i in range(len(counting_survivors)):
        if counting_survivors[best]<=counting_survivors[i]:
            best = i

    winners.append(strains[best].alpha)
    return winners


strains_id = np.linspace(0, 10, 11)
alphas = np.linspace(0, 1, 11)
for i in range(len(strains_id)):
        strain = Strain()
        print(strains)
        strain.alpha = alphas[i]
        strains.append(strain)
print('Simulation_and_Selection', Simulation_and_Selection([0,1]))


def plot_repartition(nbstrains, density='Normale'):      # parameters est un tableau car les différentes lois n'ont pas besoin du même nombre de variables
    global winners
    global strains

    #tmax : t = 2 × 108 h
    """ tmax = 2*108
    dt = 1 #h """

    new_strains = [] #Important, pour ne pas avoir à recréer les "fiches d'identités" des strains (en voulant éviter les bugs de non-effacement des tailles de populations) (oui ça pourrait être global mais jpp)

    #strains_id = np.linspace(0, 1000, 1001)
    strains_id = np.linspace(0, nbstrains, nbstrains+1)
    alphas = np.linspace(0, 1, nbstrains+1)
    for i in range(len(strains_id)):
        strain = Strain()
        strain.alpha = alphas[i] ### alphas[strains_id[i]] ?
        new_strains.append(strain)

    if density == 'Normale' :
        mus = np.linspace(0, 40, 41)*10   #de 0 à 400 heures (voir figure 4.C)
        sigmas = np.linspace(0, 100, 6)
        best_genotypes_sigma = []
        for sigma in sigmas:
            for mu in mus:
                winners = []
                parameters_proba = [mu, sigma]
                strains = copy(new_strains)
                Simulation_and_Selection(parameters_proba)
                best_genotypes_mu = winners
            best_genotypes_sigma.append(best_genotypes_mu)
        
        ordered_parameters = [] #On range les mus et sigmas selon la même structure que best_genotypes_sigma pour pouvoir tracer les courbes
        for i in range(len(sigmas)):
            subtable=[]
            for j in range(len(mus)):
                subtable.append([sigmas[i], mus[j]])
            ordered_parameters.append(subtable)
        # Il n'y a plus qu'à faire correspondre [sigma, mu==Mean_starvation] <-> [liste de genotypes]

    return best_genotypes_sigma

#print(plot_repartition(1))

""" def plot_popstrain(strain_name):
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

plot_popstrain('Strain basic') """

print("helloworld")
