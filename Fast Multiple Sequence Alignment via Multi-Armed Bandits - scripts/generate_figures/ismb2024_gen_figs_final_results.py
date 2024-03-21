# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 00:48:22 2024

@author: kmazo
"""

import numpy as np
from matplotlib import pyplot as plt

x16SBALL_spe = {'UPP': 0.052, 'UPP2': 0.043, 'UPP-fp5-K20': 0.053000000000000005}
x16SBALL_tc = {'UPP': 0.019, 'UPP2': 0.001, 'UPP-fp5-K20': 0.019}
x16SBALL_t = {'UPP': 20371.0, 'UPP2': 9139.0, 'UPP-fp5-K20': 1897.0}
x16SBALL_sps = {'UPP': 0.947, 'UPP2': 0.955, 'UPP-fp5-K20': 0.944}
x16SBALL_mod = {'UPP': 0.949, 'UPP2': 0.959, 'UPP-fp5-K20': 0.95}

x16ST_spe = {'UPP': 0.1775, 'UPP2': 0.1965, 'UPP-fp5-K20': 0.1975}
x16ST_tc = {'UPP': 0.011, 'UPP2': 0.005, 'UPP-fp5-K20': 0.009}
x16ST_t = {'UPP': 7422.0, 'UPP2': 5179.0, 'UPP-fp5-K20': 2759.0}
x16ST_sps = {'UPP': 0.831, 'UPP2': 0.792, 'UPP-fp5-K20': 0.789}
x16ST_mod = {'UPP': 0.8140000000000001, 'UPP2': 0.815, 'UPP-fp5-K20': 0.8160000000000001}



x16S3_spe = {'UPP': 0.122, 'UPP2': 0.127, 'UPP-fp5-K20': 0.1225}
x16S3_tc = {'UPP': 0.006, 'UPP2': 0.004, 'UPP-fp5-K20': 0.008}
x16S3_t = {'UPP': 6710.0, 'UPP2': 4531.0, 'UPP-fp5-K20': 2605.0}
x16S3_sps = {'UPP': 0.924, 'UPP2': 0.914, 'UPP-fp5-K20': 0.923}
x16S3_mod = {'UPP': 0.832, 'UPP2': 0.832, 'UPP-fp5-K20': 0.832}

x10000M2_spe = {'UPP': 0.0755, 'UPP2': 0.0605, 'UPP-fp5-K20': 0.062}
x10000M2_tc = {'UPP': 0.02, 'UPP2': 0.016, 'UPP-fp5-K20': 0.007}
x10000M2_t = {'UPP': 4718.0, 'UPP2': 2988.0, 'UPP-fp5-K20': 2477.0}
x10000M2_sps = {'UPP': 0.908, 'UPP2': 0.927, 'UPP-fp5-K20': 0.924}
x10000M2_mod = {'UPP': 0.9410000000000001, 'UPP2': 0.952, 'UPP-fp5-K20': 0.952}

x10000M3_spe = {'UPP': 0.0085, 'UPP2': 0.008, 'UPP-fp5-K20': 0.01}
x10000M3_tc = {'UPP': 0.113, 'UPP2': 0.077, 'UPP-fp5-K20': 0.062}
x10000M3_t = {'UPP': 3484.0, 'UPP2': 1502.0, 'UPP-fp5-K20': 3208.0}
x10000M3_sps = {'UPP': 0.988, 'UPP2': 0.988, 'UPP-fp5-K20': 0.985}
x10000M3_mod = {'UPP': 0.995, 'UPP2': 0.996, 'UPP-fp5-K20': 0.995}


x10000M4_spe = {'UPP': 0.003, 'UPP2': 0.0035, 'UPP-fp5-K20': 0.007}
x10000M4_tc = {'UPP': 0.395, 'UPP2': 0.411, 'UPP-fp5-K20': 0.115}
x10000M4_t = {'UPP': 3853.0, 'UPP2': 1470.0, 'UPP-fp5-K20': 3342.0}
x10000M4_sps =  {'UPP': 0.996, 'UPP2': 0.995, 'UPP-fp5-K20': 0.99}
x10000M4_mod = {'UPP': 0.998, 'UPP2': 0.998, 'UPP-fp5-K20': 0.996}


x10000_spe = {'UPP': 0.0955, 'UPP2': 0.0965, 'UPP-fp5-K20': 0.10600000000000001}
x10000_tc = {'UPP': 0.003, 'UPP2': 0.004, 'UPP-fp5-K20': 0.003}
x10000_t = {'UPP': 11015.0, 'UPP2': 6681.0, 'UPP-fp5-K20': 3213.0}
x10000_sps = {'UPP': 0.903, 'UPP2': 0.902, 'UPP-fp5-K20': 0.887}
x10000_mod = {'UPP': 0.906, 'UPP2': 0.905, 'UPP-fp5-K20': 0.901}


x50000_spe = {'UPP': 0.0985, 'UPP2': 0.1045, 'UPP-fp5-K20': 0.1115}
x50000_tc = {'UPP': 0.002, 'UPP2': 0.002, 'UPP-fp5-K20': 0.001}
x50000_t = {'UPP': 48182.0, 'UPP2': 23445.0, 'UPP-fp5-K20': 3986.0}
x50000_sps = {'UPP': 0.9, 'UPP2': 0.894, 'UPP-fp5-K20': 0.883}
x50000_mod = {'UPP': 0.903, 'UPP2': 0.897, 'UPP-fp5-K20': 0.894}

x100000_spe = {'UPP': 0.0895, 'UPP2': 0.0905, 'UPP-fp5-K20': 0.10850000000000001}
x100000_tc = {'UPP': 0.002, 'UPP2': 0.003, 'UPP-fp5-K20': 0.002}
x100000_t = {'UPP': 101853.0, 'UPP2': 43334.0, 'UPP-fp5-K20': 7168.0}
x100000_sps = {'UPP': 0.909, 'UPP2': 0.908, 'UPP-fp5-K20': 0.884}
x100000_mod = {'UPP': 0.912, 'UPP2': 0.911, 'UPP-fp5-K20': 0.899}


# Backbone 500


x1000M1_spe = {'UPP': 0.1895, 'UPP2': 0.4545, 'UPP-fp5-exact-K15': 0.2605, 'UPP-fp5-K20': 0.2115}

x1000M1_tc = {'UPP': 0.037, 'UPP2': 0.01, 'UPP-fp5-exact-K15': 0.005, 'UPP-fp5-K20': 0.009}

x1000M1_t = {'UPP': 1497.0, 'UPP2': 1383.0, 'UPP-fp5-exact-K15': 1454.0, 'UPP-fp5-K20': 1405.0}

x1000M1_sps = {'UPP': 0.8069999999999999, 'UPP2': 0.5389999999999999, 'UPP-fp5-exact-K15': 0.734, 'UPP-fp5-K20': 0.784}

x1000M1_mod = {'UPP': 0.8140000000000001, 'UPP2': 0.552, 'UPP-fp5-exact-K15': 0.745, 'UPP-fp5-K20': 0.793}




x1000S1_spe =  {'UPP': 0.1265, 'UPP2': 0.1905, 'UPP-fp5-exact-K15': 0.2015, 'UPP-fp5-K20': 0.17099999999999999}

x1000S1_tc = {'UPP': 0.012, 'UPP2': 0.001, 'UPP-fp5-exact-K15': 0.0, 'UPP-fp5-K20': 0.0}

x1000S1_t = {'UPP': 1289.0, 'UPP2': 1022.0, 'UPP-fp5-exact-K15': 1172.0, 'UPP-fp5-K20': 1061.0}

x1000S1_sps = {'UPP': 0.871, 'UPP2': 0.8069999999999999, 'UPP-fp5-exact-K15': 0.795, 'UPP-fp5-K20': 0.825}

x1000S1_mod = {'UPP': 0.876, 'UPP2': 0.812, 'UPP-fp5-exact-K15': 0.802, 'UPP-fp5-K20': 0.833}



x1000L1_spe = {'UPP': 0.163, 'UPP2': 0.323, 'UPP-fp5-exact-K15': 0.24, 'UPP-fp5-K20': 0.1875}

x1000L1_tc = {'UPP': 0.074, 'UPP2': 0.027, 'UPP-fp5-exact-K15': 0.017, 'UPP-fp5-K20': 0.024}

x1000L1_t = {'UPP': 1354.0, 'UPP2': 1681.0, 'UPP-fp5-exact-K15': 1334.0, 'UPP-fp5-K20': 1145.0}

x1000L1_sps = {'UPP': 0.832, 'UPP2': 0.669, 'UPP-fp5-exact-K15': 0.755, 'UPP-fp5-K20': 0.8089999999999999}

x1000L1_mod = {'UPP': 0.842, 'UPP2': 0.685, 'UPP-fp5-exact-K15': 0.765, 'UPP-fp5-K20': 0.8160000000000001}



xhomfam_spe = {'UPP': 0.24102631578947364, 'UPP2': 0.24584210526315795, 'UPP-fp5-K10': 0.25352631578947366}
xhomfam_tc =  {'UPP': 0.45952631578947367, 'UPP2': 0.4697368421052631, 'UPP-fp5-K10': 0.43768421052631573}
xhomfam_t = {'UPP': 356.7368421052632, 'UPP2': 197.1578947368421, 'UPP-fp5-K10': 272.42105263157896}
xhomfam_sps =  {'UPP': 0.873, 'UPP2': 0.923, 'UPP-fp5-K10': 0.876}
xhomfam_mod =  {'UPP': 0.780421052631579, 'UPP2': 0.7681578947368422, 'UPP-fp5-K10': 0.7808421052631579}



# x1000M1_spe = {'UPP': 0.3965, 'UPP2': 0.483, 'UPP-fp5-exact-K5': 0.8505, 'UPP-fp5-exact-K10': 0.7635000000000001, 'UPP-fp5-exact-K15': 0.3825, 'UPP-fp5-exact-K20': 0.3, 'UPP-fp5-exact-K25': 0.2695}
# x1000M1_tc = {'UPP': 0.002, 'UPP2': 0.002, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.0, 'UPP-fp5-exact-K15': 0.0, 'UPP-fp5-exact-K20': 0.0, 'UPP-fp5-exact-K25': 0.0}
# x1000M1_t = 
# x1000M1_sps = {'UPP': 0.596, 'UPP2': 0.509, 'UPP-fp5-exact-K5': 0.14500000000000002, 'UPP-fp5-exact-K10': 0.22999999999999998, 'UPP-fp5-exact-K15': 0.609, 'UPP-fp5-exact-K20': 0.6950000000000001, 'UPP-fp5-exact-K25': 0.728}
# x1000M1_mod = {'UPP': 0.611, 'UPP2': 0.525, 'UPP-fp5-exact-K5': 0.15400000000000003, 'UPP-fp5-exact-K10': 0.243, 'UPP-fp5-exact-K15': 0.626, 'UPP-fp5-exact-K20': 0.7050000000000001, 'UPP-fp5-exact-K25': 0.733}


# x1000S1_spe = {'UPP': 0.2595, 'UPP2': 0.312, 'UPP-fp5-exact-K5': 0.8194999999999999, 'UPP-fp5-exact-K10': 0.6685000000000001, 'UPP-fp5-exact-K15': 0.2945, 'UPP-fp5-exact-K20': 0.23299999999999998, 'UPP-fp5-exact-K25': 0.224}
# x1000S1_tc = {'UPP': 0.001, 'UPP2': 0.0, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.0, 'UPP-fp5-exact-K15': 0.001, 'UPP-fp5-exact-K20': 0.001, 'UPP-fp5-exact-K25': 0.0}
# x1000S1_t = 
# x1000S1_sps = {'UPP': 0.736, 'UPP2': 0.683, 'UPP-fp5-exact-K5': 0.17600000000000005, 'UPP-fp5-exact-K10': 0.32399999999999995, 'UPP-fp5-exact-K15': 0.7010000000000001, 'UPP-fp5-exact-K20': 0.763, 'UPP-fp5-exact-K25': 0.773}
# x1000S1_mod = {'UPP': 0.745, 'UPP2': 0.6930000000000001, 'UPP-fp5-exact-K5': 0.18500000000000005, 'UPP-fp5-exact-K10': 0.33899999999999997, 'UPP-fp5-exact-K15': 0.71, 'UPP-fp5-exact-K20': 0.771, 'UPP-fp5-exact-K25': 0.779}


# x1000L1_spe = {'UPP': 0.369, 'UPP2': 0.33799999999999997, 'UPP-fp5-exact-K5': 0.931, 'UPP-fp5-exact-K10': 0.7995000000000001, 'UPP-fp5-exact-K15': 0.2995, 'UPP-fp5-exact-K20': 0.261, 'UPP-fp5-exact-K25': 0.517}
# x1000L1_tc = {'UPP': 0.002, 'UPP2': 0.003, 'UPP-fp5-exact-K5': 0.0, 'UPP-fp5-exact-K10': 0.0, 'UPP-fp5-exact-K15': 0.0, 'UPP-fp5-exact-K20': 0.007, 'UPP-fp5-exact-K25': 0.003}
# x1000L1_t = 
# x1000L1_sps = {'UPP': 0.624, 'UPP2': 0.655, 'UPP-fp5-exact-K5': 0.06599999999999995, 'UPP-fp5-exact-K10': 0.19199999999999995, 'UPP-fp5-exact-K15': 0.696, 'UPP-fp5-exact-K20': 0.735, 'UPP-fp5-exact-K25': 0.469}
# x1000L1_mod = {'UPP': 0.638, 'UPP2': 0.669, 'UPP-fp5-exact-K5': 0.07199999999999995, 'UPP-fp5-exact-K10': 0.20899999999999996, 'UPP-fp5-exact-K15': 0.7050000000000001, 'UPP-fp5-exact-K20': 0.743, 'UPP-fp5-exact-K25': 0.497}

alg_names = {'UPP': 'UPP', 'UPP2': 'UPP2', 'UPP-fp5-K20': 'J-bandit'}
round_dec = 3
algorithms = x16SBALL_spe.keys()
table_str = ""
for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x16SBALL_t[alg])) +  " & " +  str(round(x16SBALL_spe[alg], round_dec)) + " & " + str(round(x16SBALL_tc[alg], round_dec)) + " & "+  str(round(x16SBALL_sps[alg], round_dec)) + " & " + str(round(x16SBALL_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x16ST_t[alg])) +  " & " +  str(round(x16ST_spe[alg], round_dec)) + " & " + str(round(x16ST_tc[alg], round_dec)) + " & "+  str(round(x16ST_sps[alg], round_dec)) + " & " + str(round(x16ST_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x16S3_t[alg])) +  " & " +  str(round(x16S3_spe[alg], round_dec)) + " & " + str(round(x16S3_tc[alg], round_dec)) + " & "+  str(round(x16S3_sps[alg], round_dec)) + " & " + str(round(x16S3_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x10000M2_t[alg])) +  " & " +  str(round(x10000M2_spe[alg], round_dec)) + " & " + str(round(x10000M2_tc[alg], round_dec)) + " & "+  str(round(x10000M2_sps[alg], round_dec)) + " & " + str(round(x10000M2_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x10000M3_t[alg])) +  " & " +  str(round(x10000M3_spe[alg], round_dec)) + " & " + str(round(x10000M3_tc[alg], round_dec)) + " & "+  str(round(x10000M3_sps[alg], round_dec)) + " & " + str(round(x10000M3_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x10000M4_t[alg])) +  " & " +  str(round(x10000M4_spe[alg], round_dec)) + " & " + str(round(x10000M4_tc[alg], round_dec)) + " & "+  str(round(x10000M4_sps[alg], round_dec)) + " & " + str(round(x10000M4_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x10000_t[alg])) +  " & " +  str(round(x10000_spe[alg], round_dec)) + " & " + str(round(x10000_tc[alg], round_dec)) + " & "+  str(round(x10000_sps[alg], round_dec)) + " & " + str(round(x10000_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x50000_t[alg])) +  " & " +  str(round(x50000_spe[alg], round_dec)) + " & " + str(round(x50000_tc[alg], round_dec)) + " & "+  str(round(x50000_sps[alg], round_dec)) + " & " + str(round(x50000_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x100000_t[alg])) +  " & " +  str(round(x100000_spe[alg], round_dec)) + " & " + str(round(x100000_tc[alg], round_dec)) + " & "+  str(round(x100000_sps[alg], round_dec)) + " & " + str(round(x100000_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x1000S1_t[alg])) +  " & " +  str(round(x1000S1_spe[alg], round_dec)) + " & " + str(round(x1000S1_tc[alg], round_dec)) + " & "+  str(round(x1000S1_sps[alg], round_dec)) + " & " + str(round(x1000S1_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x1000M1_t[alg])) +  " & " +  str(round(x1000M1_spe[alg], round_dec)) + " & " + str(round(x1000M1_tc[alg], round_dec)) + " & "+  str(round(x1000M1_sps[alg], round_dec)) + " & " + str(round(x1000M1_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

for alg in algorithms:
    table_str = table_str + " & " + alg_names[alg] + " & " + str(int(x1000L1_t[alg])) +  " & " +  str(round(x1000L1_spe[alg], round_dec)) + " & " + str(round(x1000L1_tc[alg], round_dec)) + " & "+  str(round(x1000L1_sps[alg], round_dec)) + " & " + str(round(x1000L1_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

alg_names_homfam = {'UPP': 'UPP', 'UPP2': 'UPP2', 'UPP-fp5-K10': 'J-bandit'}
algorithms_homfam = list(alg_names_homfam.keys())

for alg in algorithms_homfam:
    table_str = table_str + " & " + alg_names_homfam[alg] + " & " + str(int(xhomfam_t[alg])) +  " & " +  str(round(xhomfam_spe[alg], round_dec)) + " & " + str(round(xhomfam_tc[alg], round_dec)) + " & "+  str(round(xhomfam_sps[alg], round_dec)) + " & " + str(round(xhomfam_mod[alg], round_dec)) + " \\\ \n"
table_str = table_str + "\n"

print(table_str)



# time

alg_names1 = [alg_names[alg] for alg in alg_names]
alg_colors1 = ['red', 'blue', 'green']
fs1 = 16
fs2 = 14

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [x16SBALL_t[alg] for alg in alg_names], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("16S.B.ALL", fontsize = fs1)
plt.ylim([0, 30000])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [x16ST_t[alg] for alg in alg_names], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("16S.T", fontsize = fs1)
plt.ylim([0, 30000])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [x16S3_t[alg] for alg in alg_names], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("16S.3", fontsize = fs1)
plt.ylim([0, 30000])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [x10000M2_t[alg] for alg in alg_names], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("IND-10000M2", fontsize = fs1)
plt.ylim([0, 5000])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [x50000_t[alg] for alg in alg_names], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("RNASim-50000", fontsize = fs1)
plt.ylim([0, 110000])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [x100000_t[alg] for alg in alg_names], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("RNASim-100000", fontsize = fs1)
plt.ylim([0, 110000])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [x1000L1_t[alg] for alg in alg_names], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("ROSE-1000L1", fontsize = fs1)
plt.ylim([0, 5000])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()

fig = plt.figure( dpi=300)
plt.bar(alg_names1, [xhomfam_t[alg] for alg in alg_names_homfam], color = alg_colors1)
#plt.ylabel("time (s)")
plt.title("homfam (19)", fontsize = fs1)
plt.ylim([0, 500])
plt.xticks(fontsize = fs1)
plt.yticks(fontsize = fs2)
plt.ylabel("Time (s)", fontsize = fs1)
plt.show()




# sp score

# plt.bar(alg_names1, [x16SBALL_spe[alg] for alg in alg_names], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("16.B.ALL")
# plt.ylim([0, 1])
# plt.show()

# plt.bar(alg_names1, [x16ST_spe[alg] for alg in alg_names], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("16.T")
# plt.ylim([0, 1])
# plt.show()

# plt.bar(alg_names1, [x16S3_spe[alg] for alg in alg_names], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("16.3")
# plt.ylim([0, 1])
# plt.show()

# plt.bar(alg_names1, [x50000_spe[alg] for alg in alg_names], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("RNASim-50000")
# plt.ylim([0, 1])
# plt.show()

# plt.bar(alg_names1, [x100000_spe[alg] for alg in alg_names], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("RNASim-100000")
# plt.ylim([0, 1])
# plt.show()

# plt.bar(alg_names1, [x10000M2_spe[alg] for alg in alg_names], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("IND-10000M2")
# plt.ylim([0, 1])
# plt.show()

# plt.bar(alg_names1, [x1000L1_spe[alg] for alg in alg_names], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("ROSE-1000L1")
# plt.ylim([0, 1])
# plt.show()

# plt.bar(alg_names1, [xhomfam_spe[alg] for alg in alg_names_homfam], color = alg_colors1)
# #plt.ylabel("time (s)")
# plt.title("homfam (19)")
# plt.ylim([0, 1])
# plt.show()


# x = np.linspace(0, 2 * np.pi, 400)
# y = np.sin(x ** 2)
fig, axs = plt.subplots(2, 4, sharex=True, sharey=True, dpi=1200)
axs[0, 0].bar(alg_names1, [x16SBALL_spe[alg] for alg in alg_names], color = alg_colors1)
axs[0, 0].set_title('16S.B.ALL', fontsize =10)
axs[0, 1].bar(alg_names1, [x16ST_spe[alg] for alg in alg_names], color = alg_colors1)
axs[0, 1].set_title('16S.T', fontsize =10)
axs[0, 2].bar(alg_names1, [x16S3_spe[alg] for alg in alg_names], color = alg_colors1)
axs[0, 2].set_title('16S.3', fontsize =10)
axs[0, 3].bar(alg_names1, [x10000M2_spe[alg] for alg in alg_names], color = alg_colors1)
axs[0, 3].set_title("IND-10000M2", fontsize =10)
axs[1, 0].bar(alg_names1, [x100000_spe[alg] for alg in alg_names], color = alg_colors1)
axs[1, 0].set_title('RNASim-50000', fontsize =10)
axs[1, 1].bar(alg_names1, [x50000_spe[alg] for alg in alg_names], color = alg_colors1)
axs[1, 1].set_title('RNASim-100000', fontsize =10)
axs[1, 2].bar(alg_names1, [x1000L1_spe[alg] for alg in alg_names], color = alg_colors1)
axs[1, 2].set_title("ROSE-1000L1",  fontsize =10)
axs[1, 3].bar(alg_names1, [xhomfam_spe[alg] for alg in alg_names_homfam], color = alg_colors1)
axs[1, 3].set_title("homfam (19)", fontsize =10)
#axs[1, 0].plot(x, -y, 'tab:green')
#axs[1, 0].set_title('Axis [1, 0]')
#axs[1, 1].plot(x, -y, 'tab:red')
#axs[1, 1].set_title('Axis [1, 1]')

for ax in axs.flat:
    #ax.set(xlabel='x-label', ylabel='y-label')
    ax.set_ylabel('sp-error', fontsize = 10)
    ax.set_xticklabels(alg_names1, rotation=45, fontsize =10)
    ax.set_yticklabels([0, 0.2], fontsize =10)
    #ax.set_xticklabels(fontsize =10)

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
fig.show()