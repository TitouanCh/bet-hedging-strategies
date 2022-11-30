import numpy as np
from copy import copy
import matplotlib.pyplot as plt

class Strain:
    def __init__(self):
        self.name = "Strain basic"
        self.c = 0.173        #division rate
        self.Q = 1         #a voir, cmax/c
        self.spores = 0
        self.vg = 10**8
        self.X = self.spores + self.vg        #population
        self.alpha = 0.5
    def __str__(self):
        return "Name : " + self.name + ", Population X = " + str(self.X) + ", alpha : " + str(self.alpha)
    def __repr__(self):
        return self.__str__()
    def modifparam(self, investment):
        self.alpha = investment

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
    nu = 1.1 - 2*souche.c
    beta = 3.1 - 4*souche.c
    return (np.exp((-nu*time_starvation)**(beta)) - np.exp((-nu*T_sur)**(beta)))/(1 - np.exp((-nu*T_sur)**(beta)))

def initialisation() :
    genotype = np.linspace(0,1,1001)
    for i in range(0, 1001):
        a = Strain()
        g = genotype[i]
        a.modifparam(g)
        strains.append(a)


def starvation_step(dt):
    global Time_passed, Time_starvation, Time_rich

    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    data.append(state)

    Time_rich = 0

    for s in strains: 
        #Calcul des conséquences de l'arrivée de la starvation
        if Time_starvation == 0:            
            s.spores = s.alpha * spore_stalk * s.vg
            s.vg = (1 - s.alpha) * s.vg

        #Evolution durant starvation
        s.spores -= delta * s.spores        
        s.vg -= survival_vg(Time_starvation, s) * s.vg

    Time_starvation += dt
    Time_passed += dt


def rich_step(dt):
    global Time_passed, R, tau_counter, Time_rich, Time_starvation
    #print(Time_passed)

    Time_starvation = 0 #! ici ça doit bien être temps de starvation
    
    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    data.append(state)

    #spores germination
    if tau_counter >= tau and state['Time_passed']!=0:
        print("Germination")
        for s in strains:
            #On calcule le flux de nb de spores germinant en vg
            nu = 1.1 - 2*s.c
            s.vg += s.spores * (1 - nu)
            s.spores = 0


            ###Sens inverse :
            # flow = s.vg * s.alpha #!!! alpha bien par rapport à vg slmt ?
            # s.vg -= flow #On enlève la part correspondant aux cellules du sporocarpe, qui ne deviennent pas des endospores et ne se reproduisent pas (donc négligeables)
            # s.spores += flow * spore_stalk

        tau_counter = -1
    
    elif tau_counter >= tau : #(and Time_passed == 0)
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
    
    Time_rich += dt




#fig4

# B
#def starvtime_exp():

    t = np.linspace(0, 400,401)
    plt.plot(t, ((1/50)*np.exp(-(1/50)*t)))
    plt.show()




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
