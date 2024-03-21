# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 22:33:58 2024

@author: kmazo
"""
import numpy as np
from matplotlib import pyplot as plt

x16SBALL_spe = {'UPP': 0.052, 'UPP2': 0.043, 'UPP-fp5-K5': 0.0645, 'UPP-fp5-K10': 0.0535, 'UPP-fp5-K15': 0.05, 'UPP-fp5-K20': 0.051500000000000004, 'UPP-fp5-K25': 0.052500000000000005, 'UPP-fp5-exact-K5': 0.068, 'UPP-fp5-exact-K10': 0.05, 'UPP-fp5-exact-K15': 0.051500000000000004, 'UPP-fp5-exact-K20': 0.0495, 'UPP-fp5-exact-K25': 0.048}
x16SBALL_tc = {'UPP': 0.019, 'UPP2': 0.001, 'UPP-fp5-K5': 0.01, 'UPP-fp5-K10': 0.018, 'UPP-fp5-K15': 0.019, 'UPP-fp5-K20': 0.023, 'UPP-fp5-K25': 0.018, 'UPP-fp5-exact-K5': 0.02, 'UPP-fp5-exact-K10': 0.018, 'UPP-fp5-exact-K15': 0.018, 'UPP-fp5-exact-K20': 0.017, 'UPP-fp5-exact-K25': 0.019}
x16SBALL_sps = {'UPP': 0.947, 'UPP2': 0.955, 'UPP-fp5-K5': 0.931, 'UPP-fp5-K10': 0.944, 'UPP-fp5-K15': 0.948, 'UPP-fp5-K20': 0.945, 'UPP-fp5-K25': 0.943, 'UPP-fp5-exact-K5': 0.926, 'UPP-fp5-exact-K10': 0.948, 'UPP-fp5-exact-K15': 0.947, 'UPP-fp5-exact-K20': 0.948, 'UPP-fp5-exact-K25': 0.951}
x16SBALL_mod = {'UPP': 0.949, 'UPP2': 0.959, 'UPP-fp5-K5': 0.94, 'UPP-fp5-K10': 0.949, 'UPP-fp5-K15': 0.952, 'UPP-fp5-K20': 0.952, 'UPP-fp5-K25': 0.952, 'UPP-fp5-exact-K5': 0.938, 'UPP-fp5-exact-K10': 0.952, 'UPP-fp5-exact-K15': 0.95, 'UPP-fp5-exact-K20': 0.953, 'UPP-fp5-exact-K25': 0.953}


x10000M2_spe = {'UPP': 0.0755, 'UPP2': 0.0605, 'UPP-fp5-exact-K5': 0.953, 'UPP-fp5-exact-K10': 0.7909999999999999, 'UPP-fp5-exact-K15': 0.5745, 'UPP-fp5-exact-K20': 0.057999999999999996, 'UPP-fp5-exact-K25': 0.0565}
x10000M2_tc = {'UPP': 0.02, 'UPP2': 0.016, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.002, 'UPP-fp5-exact-K15': 0.0, 'UPP-fp5-exact-K20': 0.021, 'UPP-fp5-exact-K25': 0.017}
x10000M2_sps = {'UPP': 0.908, 'UPP2': 0.927, 'UPP-fp5-exact-K5': 0.03200000000000003, 'UPP-fp5-exact-K10': 0.11499999999999999, 'UPP-fp5-exact-K15': 0.254, 'UPP-fp5-exact-K20': 0.927, 'UPP-fp5-exact-K25': 0.931}
x10000M2_mod = {'UPP': 0.9410000000000001, 'UPP2': 0.952, 'UPP-fp5-exact-K5': 0.062000000000000055, 'UPP-fp5-exact-K10': 0.30300000000000005, 'UPP-fp5-exact-K15': 0.597, 'UPP-fp5-exact-K20': 0.957, 'UPP-fp5-exact-K25': 0.956}


x10000M3_spe = {'UPP': 0.0085, 'UPP2': 0.008, 'UPP-fp5-exact-K5': 0.609, 'UPP-fp5-exact-K10': 0.4345, 'UPP-fp5-exact-K15': 0.35150000000000003, 'UPP-fp5-exact-K20': 0.0105, 'UPP-fp5-exact-K25': 0.009}
x10000M3_tc = {'UPP': 0.113, 'UPP2': 0.077, 'UPP-fp5-exact-K5': 0.001, 'UPP-fp5-exact-K10': 0.024, 'UPP-fp5-exact-K15': 0.008, 'UPP-fp5-exact-K20': 0.069, 'UPP-fp5-exact-K25': 0.077}
x10000M3_sps = {'UPP': 0.988, 'UPP2': 0.988, 'UPP-fp5-exact-K5': 0.271, 'UPP-fp5-exact-K10': 0.395, 'UPP-fp5-exact-K15': 0.486, 'UPP-fp5-exact-K20': 0.984, 'UPP-fp5-exact-K25': 0.987}
x10000M3_mod = {'UPP': 0.995, 'UPP2': 0.996, 'UPP-fp5-exact-K5': 0.511, 'UPP-fp5-exact-K10': 0.736, 'UPP-fp5-exact-K15': 0.8109999999999999, 'UPP-fp5-exact-K20': 0.995, 'UPP-fp5-exact-K25': 0.995}


x10000M4_spe = {'UPP': 0.003, 'UPP2': 0.0035, 'UPP-fp5-exact-K5': 0.0715, 'UPP-fp5-exact-K10': 0.02, 'UPP-fp5-exact-K15': 0.0305, 'UPP-fp5-exact-K20': 0.01, 'UPP-fp5-exact-K25': 0.003}
x10000M4_tc = {'UPP': 0.395, 'UPP2': 0.411, 'UPP-fp5-exact-K5': 0.027, 'UPP-fp5-exact-K10': 0.121, 'UPP-fp5-exact-K15': 0.072, 'UPP-fp5-exact-K20': 0.101, 'UPP-fp5-exact-K25': 0.277}
x10000M4_sps =  {'UPP': 0.996, 'UPP2': 0.995, 'UPP-fp5-exact-K5': 0.911, 'UPP-fp5-exact-K10': 0.974, 'UPP-fp5-exact-K15': 0.961, 'UPP-fp5-exact-K20': 0.987, 'UPP-fp5-exact-K25': 0.995}
x10000M4_mod = {'UPP': 0.998, 'UPP2': 0.998, 'UPP-fp5-exact-K5': 0.946, 'UPP-fp5-exact-K10': 0.986, 'UPP-fp5-exact-K15': 0.978, 'UPP-fp5-exact-K20': 0.993, 'UPP-fp5-exact-K25': 0.999}


x10000_spe = {'UPP': 0.0955, 'UPP2': 0.0965, 'UPP-fp5-exact-K5': 0.1855, 'UPP-fp5-exact-K10': 0.1065, 'UPP-fp5-exact-K15': 0.1085, 'UPP-fp5-exact-K20': 0.10600000000000001, 'UPP-fp5-exact-K25': 0.1205}
x10000_tc = {'UPP': 0.003, 'UPP2': 0.004, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.006, 'UPP-fp5-exact-K15': 0.004, 'UPP-fp5-exact-K20': 0.003, 'UPP-fp5-exact-K25': 0.004}
x10000_sps = {'UPP': 0.903, 'UPP2': 0.902, 'UPP-fp5-exact-K5': 0.788, 'UPP-fp5-exact-K10': 0.891, 'UPP-fp5-exact-K15': 0.888, 'UPP-fp5-exact-K20': 0.887, 'UPP-fp5-exact-K25': 0.866}
x10000_mod = {'UPP': 0.906, 'UPP2': 0.905, 'UPP-fp5-exact-K5': 0.841, 'UPP-fp5-exact-K10': 0.896, 'UPP-fp5-exact-K15': 0.895, 'UPP-fp5-exact-K20': 0.901, 'UPP-fp5-exact-K25': 0.893}


x1000M1_spe = {'UPP': 0.3965, 'UPP2': 0.483, 'UPP-fp5-exact-K5': 0.8505, 'UPP-fp5-exact-K10': 0.7635000000000001, 'UPP-fp5-exact-K15': 0.3825, 'UPP-fp5-exact-K20': 0.3, 'UPP-fp5-exact-K25': 0.2695}
x1000M1_tc = {'UPP': 0.002, 'UPP2': 0.002, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.0, 'UPP-fp5-exact-K15': 0.0, 'UPP-fp5-exact-K20': 0.0, 'UPP-fp5-exact-K25': 0.0}
x1000M1_sps = {'UPP': 0.596, 'UPP2': 0.509, 'UPP-fp5-exact-K5': 0.14500000000000002, 'UPP-fp5-exact-K10': 0.22999999999999998, 'UPP-fp5-exact-K15': 0.609, 'UPP-fp5-exact-K20': 0.6950000000000001, 'UPP-fp5-exact-K25': 0.728}
x1000M1_mod = {'UPP': 0.611, 'UPP2': 0.525, 'UPP-fp5-exact-K5': 0.15400000000000003, 'UPP-fp5-exact-K10': 0.243, 'UPP-fp5-exact-K15': 0.626, 'UPP-fp5-exact-K20': 0.7050000000000001, 'UPP-fp5-exact-K25': 0.733}


x1000S1_spe = {'UPP': 0.2595, 'UPP2': 0.312, 'UPP-fp5-exact-K5': 0.8194999999999999, 'UPP-fp5-exact-K10': 0.6685000000000001, 'UPP-fp5-exact-K15': 0.2945, 'UPP-fp5-exact-K20': 0.23299999999999998, 'UPP-fp5-exact-K25': 0.224}
x1000S1_tc = {'UPP': 0.001, 'UPP2': 0.0, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.0, 'UPP-fp5-exact-K15': 0.001, 'UPP-fp5-exact-K20': 0.001, 'UPP-fp5-exact-K25': 0.0}
x1000S1_sps = {'UPP': 0.736, 'UPP2': 0.683, 'UPP-fp5-exact-K5': 0.17600000000000005, 'UPP-fp5-exact-K10': 0.32399999999999995, 'UPP-fp5-exact-K15': 0.7010000000000001, 'UPP-fp5-exact-K20': 0.763, 'UPP-fp5-exact-K25': 0.773}
x1000S1_mod = {'UPP': 0.745, 'UPP2': 0.6930000000000001, 'UPP-fp5-exact-K5': 0.18500000000000005, 'UPP-fp5-exact-K10': 0.33899999999999997, 'UPP-fp5-exact-K15': 0.71, 'UPP-fp5-exact-K20': 0.771, 'UPP-fp5-exact-K25': 0.779}


x1000L1_spe = {'UPP': 0.369, 'UPP2': 0.33799999999999997, 'UPP-fp5-exact-K5': 0.931, 'UPP-fp5-exact-K10': 0.7995000000000001, 'UPP-fp5-exact-K15': 0.2995, 'UPP-fp5-exact-K20': 0.261, 'UPP-fp5-exact-K25': 0.517}
x1000L1_tc = {'UPP': 0.002, 'UPP2': 0.003, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.0, 'UPP-fp5-exact-K15': 0.0, 'UPP-fp5-exact-K20': 0.007, 'UPP-fp5-exact-K25': 0.003}
x1000L1_sps = {'UPP': 0.624, 'UPP2': 0.655, 'UPP-fp5-exact-K5': 0.06599999999999995, 'UPP-fp5-exact-K10': 0.19199999999999995, 'UPP-fp5-exact-K15': 0.696, 'UPP-fp5-exact-K20': 0.735, 'UPP-fp5-exact-K25': 0.469}
x1000L1_mod = {'UPP': 0.638, 'UPP2': 0.669, 'UPP-fp5-exact-K5': 0.07199999999999995, 'UPP-fp5-exact-K10': 0.20899999999999996, 'UPP-fp5-exact-K15': 0.7050000000000001, 'UPP-fp5-exact-K20': 0.743, 'UPP-fp5-exact-K25': 0.497}



# fig = plt.figure()
# ax = plt.subplot(111)

# for i in xrange(5):
#     ax.plot(x, i * x, label='$y = %ix$'%i)

# # Shrink current axis by 20%
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# # Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))



# Ks = [5, 10, 15, 20, 25]
# K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
# plt.plot(Ks, [x16SBALL_spe[K_to_run[K]] for K in Ks])
# # plt.plot(Ks, len(Ks)*[x16SBALL_spe['UPP']])
# # plt.plot(Ks, [x16SBALL_spe['UPP']] + (len(Ks)-1)*[None])

# sp-error

colors1 = ['red', 'blue', 'green', 'purple', 'brown', 'black', 'orange', 'cyan']
markerstyle1 = "+"
markersize1 = 15
fig = plt.figure( dpi=300) #1200)
ax = plt.subplot(111)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x16SBALL_spe[K_to_run[K]] for K in Ks], label = "16S.B.ALL", color = colors1[0])
#ax.plot(Ks, [x16SBALL_spe[K_to_run[K]] for K in Ks], label = "16S.B.ALL")
plt.plot([4], [x16SBALL_spe['UPP']], marker = markerstyle1, color = colors1[0], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M2_spe[K_to_run[K]] for K in Ks], label = "10000M2", color = colors1[1])
plt.plot([4], [x10000M2_spe['UPP']], marker = markerstyle1, color = colors1[1], markersize = markersize1)


Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M3_spe[K_to_run[K]] for K in Ks], label = "10000M3", color = colors1[2])
plt.plot([4], [x10000M3_spe['UPP']], marker = markerstyle1, color = colors1[2], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M4_spe[K_to_run[K]] for K in Ks], label = "10000M4", color = colors1[3])
plt.plot([4], [x10000M4_spe['UPP']], marker = markerstyle1, color = colors1[3], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000_spe[K_to_run[K]] for K in Ks], label = "RNASim10000", color = colors1[4])
plt.plot([4], [x10000_spe['UPP']], marker = markerstyle1, color = colors1[4], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000S1_spe[K_to_run[K]] for K in Ks], label = "1000S1", color = colors1[5])
plt.plot([4], [x1000S1_spe['UPP']], marker = markerstyle1, color = colors1[5], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000M1_spe[K_to_run[K]] for K in Ks], label = "1000M1", color = colors1[6])
plt.plot([4], [x1000M1_spe['UPP']], marker = markerstyle1, color = colors1[6], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000L1_spe[K_to_run[K]] for K in Ks], label = "1000L1", color = colors1[7])
plt.plot([4], [x1000L1_spe['UPP']], marker = markerstyle1, color = colors1[7], markersize = markersize1)


box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#plt.legend()
plt.title("SP-Error versus k-mer Length")
plt.xlabel("k-mer Length")
plt.xticks([5, 10, 15, 20, 25])
plt.ylabel("SP-error")
#plt.ylim([0, 1])
plt.xlim([3, 25])
plt.show()



# total column

# colors1 = ['red', 'blue', 'green', 'purple', 'brown', 'black', 'orange', 'cyan']
# markerstyle1 = "+"
# markersize1 = 10
fig = plt.figure( dpi=300) #1200)
ax = plt.subplot(111)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x16SBALL_tc[K_to_run[K]] for K in Ks], label = "16S.B.ALL", color = colors1[0])
#ax.plot(Ks, [x16SBALL_spe[K_to_run[K]] for K in Ks], label = "16S.B.ALL")
plt.plot([4], [x16SBALL_tc['UPP']], marker = markerstyle1, color = colors1[0], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M2_tc[K_to_run[K]] for K in Ks], label = "10000M2", color = colors1[1])
plt.plot([4], [x10000M2_tc['UPP']], marker = markerstyle1, color = colors1[1], markersize = markersize1)


Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M3_tc[K_to_run[K]] for K in Ks], label = "10000M3", color = colors1[2])
plt.plot([4], [x10000M3_tc['UPP']], marker = markerstyle1, color = colors1[2], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M4_tc[K_to_run[K]] for K in Ks], label = "10000M4", color = colors1[3])
plt.plot([4], [x10000M4_tc['UPP']], marker = markerstyle1, color = colors1[3], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000_tc[K_to_run[K]] for K in Ks], label = "RNASim10000", color = colors1[4])
plt.plot([4], [x10000_tc['UPP']], marker = markerstyle1, color = colors1[4], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000S1_tc[K_to_run[K]] for K in Ks], label = "1000S1", color = colors1[5])
plt.plot([4], [x1000S1_tc['UPP']], marker = markerstyle1, color = colors1[5], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000M1_tc[K_to_run[K]] for K in Ks], label = "1000M1", color = colors1[6])
plt.plot([4], [x1000M1_tc['UPP']], marker = markerstyle1, color = colors1[6], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000L1_tc[K_to_run[K]] for K in Ks], label = "1000L1", color = colors1[7])
plt.plot([4], [x1000L1_tc['UPP']], marker = markerstyle1, color = colors1[7], markersize = markersize1)


box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#plt.legend()
plt.title("TC-Score vesus k-mer Length")
plt.xlabel("k-mer Length")
plt.xticks([5, 10, 15, 20, 25])
plt.ylabel("TC-Score")
#plt.ylim([0, 1])
plt.xlim([3, 25])
plt.show()


# sp score

# colors1 = ['red', 'blue', 'green', 'purple', 'brown', 'black', 'orange', 'cyan']
# markerstyle1 = "+"
# markersize1 = 10
fig = plt.figure( dpi=300) #1200)
ax = plt.subplot(111)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x16SBALL_sps[K_to_run[K]] for K in Ks], label = "16S.B.ALL", color = colors1[0])
#ax.plot(Ks, [x16SBALL_spe[K_to_run[K]] for K in Ks], label = "16S.B.ALL")
plt.plot([4], [x16SBALL_sps['UPP']], marker = markerstyle1, color = colors1[0], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M2_sps[K_to_run[K]] for K in Ks], label = "10000M2", color = colors1[1])
plt.plot([4], [x10000M2_sps['UPP']], marker = markerstyle1, color = colors1[1], markersize = markersize1)


Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M3_sps[K_to_run[K]] for K in Ks], label = "10000M3", color = colors1[2])
plt.plot([4], [x10000M3_sps['UPP']], marker = markerstyle1, color = colors1[2], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M4_sps[K_to_run[K]] for K in Ks], label = "10000M4", color = colors1[3])
plt.plot([4], [x10000M4_sps['UPP']], marker = markerstyle1, color = colors1[3], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000_sps[K_to_run[K]] for K in Ks], label = "RNASim10000", color = colors1[4])
plt.plot([4], [x10000_sps['UPP']], marker = markerstyle1, color = colors1[4], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000S1_sps[K_to_run[K]] for K in Ks], label = "1000S1", color = colors1[5])
plt.plot([4], [x1000S1_sps['UPP']], marker = markerstyle1, color = colors1[5], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000M1_sps[K_to_run[K]] for K in Ks], label = "1000M1", color = colors1[6])
plt.plot([4], [x1000M1_sps['UPP']], marker = markerstyle1, color = colors1[6], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000L1_sps[K_to_run[K]] for K in Ks], label = "1000L1", color = colors1[7])
plt.plot([4], [x1000L1_sps['UPP']], marker = markerstyle1, color = colors1[7], markersize = markersize1)


box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#plt.legend()
plt.title("SP-Score versus k-mer Length")
plt.xlabel("k-mer Length")
plt.xticks([5, 10, 15, 20, 25])
plt.ylabel("SP-Score")
#plt.ylim([0, 1])
plt.xlim([3, 25])
plt.show()



# modeler score

# colors1 = ['red', 'blue', 'green', 'purple', 'brown', 'black', 'orange', 'cyan']
# markerstyle1 = "+"
# markersize1 = 15
fig = plt.figure( dpi=300) #1200)
ax = plt.subplot(111)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x16SBALL_mod[K_to_run[K]] for K in Ks], label = "16S.B.ALL", color = colors1[0])
#ax.plot(Ks, [x16SBALL_spe[K_to_run[K]] for K in Ks], label = "16S.B.ALL")
plt.plot([4], [x16SBALL_mod['UPP']], marker = markerstyle1, color = colors1[0], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M2_mod[K_to_run[K]] for K in Ks], label = "10000M2", color = colors1[1])
plt.plot([4], [x10000M2_mod['UPP']], marker = markerstyle1, color = colors1[1], markersize = markersize1)


Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M3_mod[K_to_run[K]] for K in Ks], label = "10000M3", color = colors1[2])
plt.plot([4], [x10000M3_mod['UPP']], marker = markerstyle1, color = colors1[2], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000M4_mod[K_to_run[K]] for K in Ks], label = "10000M4", color = colors1[3])
plt.plot([4], [x10000M4_mod['UPP']], marker = markerstyle1, color = colors1[3], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x10000_mod[K_to_run[K]] for K in Ks], label = "RNASim10000", color = colors1[4])
plt.plot([4], [x10000_mod['UPP']], marker = markerstyle1, color = colors1[4], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000S1_mod[K_to_run[K]] for K in Ks], label = "1000S1", color = colors1[5])
plt.plot([4], [x1000S1_mod['UPP']], marker = markerstyle1, color = colors1[5], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000M1_mod[K_to_run[K]] for K in Ks], label = "1000M1", color = colors1[6])
plt.plot([4], [x1000M1_mod['UPP']], marker = markerstyle1, color = colors1[6], markersize = markersize1)

Ks = [5, 10, 15, 20, 25]
K_to_run = {5: 'UPP-fp5-exact-K5', 10: 'UPP-fp5-exact-K10', 15: 'UPP-fp5-exact-K15', 20: 'UPP-fp5-exact-K20', 25: 'UPP-fp5-exact-K25'}
plt.plot(Ks, [x1000L1_mod[K_to_run[K]] for K in Ks], label = "1000L1", color = colors1[7])
plt.plot([4], [x1000L1_mod['UPP']], marker = markerstyle1, color = colors1[7], markersize = markersize1)


box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

#plt.legend()
plt.title("Modeler Score versus k-mer Length")
plt.xlabel("k-mer Length")
plt.xticks([5, 10, 15, 20, 25])
plt.ylabel("Modeler Score")
#plt.ylim([0, 1])
plt.xlim([3, 25])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.show()

