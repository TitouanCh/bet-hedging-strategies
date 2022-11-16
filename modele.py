import numpy as np
from copy import copy

class Strain:
    def __init__(self):
        self.name = "Strain basic"
        self.c = 0.173        #division rate
        self.Q = 1         #a voir, cmax/c
        self.X = 1000        #population
    def __str__(self):
        return "Name : " + self.name + ", Population X = " + str(self.X)
    def __repr__(self):
        return self.__str__()

def copyStrains(strains):
    resultat = []
    for s in strains:
        resultat.append(copy(s))
    return resultat

data = []
R0 = 10**8    #ressources
Rhalf = 0.1*R0
R = R0
strains = []
Time_passed = 0       

def initialisation() :
    a = Strain()
    print(a)
    strains.append(a)

def rich_step(dt):
    global Time_passed, R
    #print(Time_passed)
    
    #stockage
    state = {"Time_passed" : Time_passed, "R" : R, "strains" : copyStrains(strains)}
    data.append(state)

    #calcul
    kRGab = state["R"] / (state["R"] + Rhalf)
    sommeTit = 0

    #calcul - population
    for s in strains:
        sommeTit += s.Q * s.c * s.X
        s.X += s.c * kRGab * s.X * dt

    # calcul - ressource
    R += - kRGab * sommeTit

    Time_passed += dt

initialisation()
for i in range(1,6) :
    rich_step(1)

print(data)

print("helloworld")
