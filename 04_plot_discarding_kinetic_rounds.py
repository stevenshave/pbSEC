import numpy as np
from high_accuracy_binding_equations import *
import matplotlib.pyplot as plt

XAXIS_BEGINNING = 0  # pKD of 3 is mM
XAXIS_END = 1  # pKD of 12 is pM
YAXIS_BEGINNING = 1e-3
YAXIS_END = 10
NUM_POINTS_ON_XAXIS = 1000  # Publication used 2000 pts along X
PROTEIN_CONC = 10.0
LIGAND_CONC = 8.0  # Singular compound conc in pool
KD = 10  # 1 µM
NUM_ROUNDS = 5

print("Hello")
x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
#ligand_kd_range = 10**(-x_axis)*1e6
ligand_concs = np.full((NUM_ROUNDS, NUM_POINTS_ON_XAXIS), np.nan)
one_round_lig_conc = one_to_one_binding(PROTEIN_CONC, LIGAND_CONC, KD)

for efficiency_i, efficiency in enumerate(x_axis):
    ligand_concs[0, efficiency_i] = one_round_lig_conc*efficiency
for round_num in range(1, NUM_ROUNDS):
    for efficiency_i, efficiency in enumerate(x_axis):
        ligand_concs[round_num, efficiency_i] = one_to_one_binding(
            PROTEIN_CONC*efficiency**round_num, ligand_concs[round_num-1, efficiency_i]*efficiency**round_num, KD)

ligand_concs*=1e-6

fig, ax = plt.subplots(1, 1, figsize=(7.204724, 5.09424929292))

for i in range(NUM_ROUNDS):
    ax.plot(x_axis, ligand_concs[i], label="Round "+str(i+1))
ax.hlines(1e-8,0,1)
plt.legend()
plt.yscale('log')
# ax.set_xticklabels(
#     ["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)"])
# ax.set_xticks(range(XAXIS_BEGINNING, XAXIS_END+1))
ax.set_xlim(0,1)
ax.set_ylim(1e-9, 1e-5)
ax.set_xlabel(r"Recovery efficiency", fontsize=16)
ax.set_ylabel(r"[Ligand] (M)", fontsize=16)
ax.set_title(r'pbSEC kinetic scheme, [Prot] = 10 $\mathrm{\mu}$M, [Ligand] = 8 $\mathrm{\mu}$M, K$_\mathrm{D}$ = '+str(KD)+' µM')
plt.show()
