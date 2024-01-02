"""
Plot pbSEC rounds, pKd vs concentration

"""

import numpy as np
from pbsec import pbSEC_simulate_n_rounds
import matplotlib.pyplot as plt

XAXIS_BEGINNING = 3  # pKD of 3 is 1 mM
XAXIS_END = 9  # pKD of 12 is 1 pM
NUM_POINTS_ON_XAXIS = 2000  # Publication used 2000 pts along X
PROTEIN_CONC = 10.0
LIGAND_CONC = 8.0  # Singular compound conc in pool
NUM_ROUNDS = 10
RECOVERY_EFFICIENCY = 1.0

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
ligand_kd_range = 10 ** (-x_axis) * 1e6
protein_concs = np.full((NUM_ROUNDS, NUM_POINTS_ON_XAXIS), np.nan)

for kd_i, kd in enumerate(ligand_kd_range):
    protein_concs[:, kd_i] = pbSEC_simulate_n_rounds(
        PROTEIN_CONC,
        LIGAND_CONC,
        kd,
        NUM_ROUNDS,
        recovery_efficiency=RECOVERY_EFFICIENCY,
    )

fig, ax = plt.subplots(1, 1, figsize=(7.204724, 5.09424929292))

for i in range(NUM_ROUNDS):
    ax.plot(x_axis, protein_concs[i] * 1e-6, label="Round " + str(i + 1))
plt.legend()
plt.yscale("log")
ax.set_xticklabels(["3 (mM)", "4", "5", r"6 ($\mathrm{\mu}$M)", "7", "8", "9 (nM)"])
ax.set_xlim(3, 9)
ax.set_xlabel(r"Ligand pK$_\mathrm{D}$", fontsize=16)
ax.set_ylabel(r"[Ligand] ($\mathrm{\mu}$M)", fontsize=16)
ax.set_title(
    "pbSEC kinetic scheme, [Prot] = 10 $\mathrm{\mu}$M, [Ligand] = 8 $\mathrm{\mu}$M"
)
plt.show()
