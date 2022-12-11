import json
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    f = open("./outputs/" + filename, "r")
    inside = json.loads(f.read())
    f.close()
    return inside

determinist = read_file("determinist.txt")
exponential = read_file("exponential.txt")
normal_20 = read_file("normal_20.txt")
normal_40 = read_file("normal_40.txt")
normal_60 = read_file("normal_60.txt")

f, (ax1, ax2) = plt.subplots(2, 2, sharex=True, figsize=(16, 9), gridspec_kw={'width_ratios': [2, 1]})

ax1[0].plot(determinist[0], determinist[1])
ax1[0].plot(exponential[0], exponential[1])

ax2[0].plot(determinist[0], determinist[1])
ax2[0].plot(normal_20[0], normal_20[1])
ax2[0].plot(normal_40[0], normal_40[1])
ax2[0].plot(normal_60[0], normal_60[1])

t = np.linspace(0, 400,401)
ax1[1].plot(t, ((1/50)*np.exp(-(1/50)*t)))

t = np.ones(400) * 50
ax2[1].plot(t, t * np.linspace(0, 0.03, 400))

normal_density = lambda x, o, mu : 1/(o*np.sqrt(2 * np.pi)) * np.exp(-1/2*np.power((x-mu)/o, 2))

t = np.linspace(0, 200, 201)
ax2[1].plot(t, normal_density(t, 20, 50))
ax2[1].plot(t, normal_density(t, 40, 50))
ax2[1].plot(t, normal_density(t, 60, 50))

ax2[1].set_ylim([0, 0.03])

plt.show()