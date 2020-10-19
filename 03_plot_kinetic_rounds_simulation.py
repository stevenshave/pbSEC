import numpy as np
from high_accuracy_binding_equations import *
import matplotlib.pyplot as plt

XAXIS_BEGINNING = 3  # pKD of 3 is mM
XAXIS_END = 9  # pKD of 12 is pM
NUM_POINTS_ON_XAXIS = 100 # Publication used 2000 pts along X
PROTEIN_CONC=10.0
LIGAND_CONC=8.0 # Singular compound conc in pool
NUM_ROUNDS=3
VOLUME=75.0 # 75 µl

RECOVERY_EFFICIENCY=0.9


print("Hello")
x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
ligand_kd_range = 10**(-x_axis)*1e6
protein_concs = np.full((NUM_ROUNDS, NUM_POINTS_ON_XAXIS) , np.nan)

for kd_i,kd in enumerate(ligand_kd_range):
    protein_concs[0,kd_i]=one_to_one_binding(PROTEIN_CONC,LIGAND_CONC,kd)
for round_num in range(1,NUM_ROUNDS):
    for kd_i,kd in enumerate(ligand_kd_range):
        protein_concs[round_num,kd_i]=one_to_one_binding(PROTEIN_CONC*RECOVERY_EFFICIENCY**round_num,protein_concs[round_num-1,kd_i]*RECOVERY_EFFICIENCY**round_num,kd)
    

fig, ax = plt.subplots(1,1, figsize=(7.204724, 5.09424929292))

for i in range(NUM_ROUNDS):
    ax.plot(x_axis, protein_concs[i], label="Round "+str(i+1))
#ax.hlines(1,3,9)
plt.legend()

ax.set_xticklabels(
    ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)"])
ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))
ax.set_xlim(3,9)
ax.set_ylim(0,8)
ax.set_xlabel(r"Ligand pK$_\mathrm{D}$", fontsize=16)
ax.set_ylabel(r"[Ligand] ($\mathrm{\mu}$M)", fontsize=16)
ax.set_title('pbSEC kinetic scheme, [Prot] = 10 $\mathrm{\mu}$M, [Ligand] = 8 $\mathrm{\mu}$M')
plt.show()

