import json
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    f = open("./outputs/" + filename, "r")
    inside = json.loads(f.read())
    f.close()
    return inside

determinist = []
for i in ["0", "3", "6", "9", "11", "12"]:
    determinist.append(read_file("severedry_d" + i + "_10-4.txt"))

exponential = []
for i in ["0", "3", "6", "9", "11", "12"]:
    exponential.append(read_file("severedry_e" + i + "_10-4.txt"))

normal_20 = []
for i in ["0", "3", "6", "9", "11", "12"]:
    normal_20.append(read_file("severedry_n20_" + i + "_10-4.txt"))

normal_60 = []
for i in ["0", "3", "6", "9", "11", "12"]:
    normal_60.append(read_file("severedry_n60_" + i + "_10-4.txt"))

f, (ax1, ax2) = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 9))

colors = ["orange", "lightgreen", "purple", "blue", "red", "darkgreen"]
label = ["$T_{dry} =$ 0", "$T_{dry} =$ 3 mois", "$T_{dry} =$ 6 mois", "$T_{dry} =$ 9 mois", "$T_{dry} =$ 11 mois", "$T_{dry} =$ 12 mois"]

for i in range(len(determinist)):
    ax1[0].plot(determinist[i][0], determinist[i][1], label=label[i], c=colors[i], linestyle='--', marker='s', markersize=3.0)

ax1[0].legend(loc="lower right")
f.supylabel('Investissement dans les spores, α')
f.supxlabel('Temps moyen de famine en saison humide $λ_{T}$ (heure)')

for i in range(len(exponential)):
    ax1[1].plot(exponential[i][0], exponential[i][1], c=colors[i], linestyle='--', marker='s', markersize=3.0)

for i in range(len(normal_20)):
    ax2[0].plot(normal_20[i][0], normal_20[i][1], c=colors[i], linestyle='--', marker='s', markersize=3.0)

for i in range(len(normal_60)):
    ax2[1].plot(normal_60[i][0], normal_60[i][1], c=colors[i], linestyle='--', marker='s', markersize=3.0)

"""
ax2[0].plot(determinist[0], determinist[1], label="σ = 0 (det.)", c="red", linestyle='--', marker='s', markersize=3.0)
ax2[0].plot(normal_20[0], normal_20[1], label="σ = 20", c="green", linestyle='--', marker='s', markersize=3.0)
ax2[0].plot(normal_40[0], normal_40[1], label="σ = 40", c="blue", linestyle='--', marker='s', markersize=3.0)
ax2[0].plot(normal_60[0], normal_60[1], label="σ = 60", c="orange", linestyle='--', marker='s', markersize=3.0)
ax2[0].legend(loc="lower right")
ax2[0].set_xlabel("Mean starvation time $λ_{T}$")
ax2[1].set_xlabel("Starvation time")

t = np.linspace(0, 400,401)
ax1[1].plot(t, ((1/50)*np.exp(-(1/50)*t)), c="black")

t = np.ones(400) * 50
ax2[1].plot(t, t * np.linspace(0, 0.03, 400), c="red")

normal_density = lambda x, o, mu : 1/(o*np.sqrt(2 * np.pi)) * np.exp(-1/2*np.power((x-mu)/o, 2))

t = np.linspace(0, 200, 201)
ax2[1].plot(t, normal_density(t, 20, 50), c="green")
ax2[1].plot(t, normal_density(t, 40, 50), c="blue")
ax2[1].plot(t, normal_density(t, 60, 50), c="orange")

ax2[1].set_ylim([0, 0.03])
"""
plt.show()