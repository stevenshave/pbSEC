"""
Plot pbSEC rounds, recovery efficiency vs log complex

"""

import numpy as np
from pbsec import pbSEC_simulate_n_rounds
import matplotlib.pyplot as plt

XAXIS_BEGINNING = 0  # pKD of 3 is mM
XAXIS_END = 100  # pKD of 12 is pM
YAXIS_BEGINNING = 1e-3
YAXIS_END = 10
NUM_POINTS_ON_XAXIS = 2000  # Publication used 2000 pts along X
PROTEIN_CONC = 10.0
LIGAND_CONC = 8.0  # Singular compound conc in pool
KD = 10  # 1 µM
NUM_ROUNDS = 5

x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
# ligand_kd_range = 10**(-x_axis)*1e6
ligand_concs = np.full((NUM_ROUNDS, NUM_POINTS_ON_XAXIS), np.nan)

for efficiency_i, efficiency in enumerate(x_axis / 100.0):
    ligand_concs[:, efficiency_i] = pbSEC_simulate_n_rounds(
        PROTEIN_CONC, LIGAND_CONC, KD, NUM_ROUNDS, recovery_efficiency=efficiency
    )

ligand_concs *= 1e-6

fig, ax = plt.subplots(1, 1, figsize=(7.204724, 5.09424929292))

for i in range(NUM_ROUNDS):
    ax.plot(x_axis, ligand_concs[i], label="Round " + str(i + 1))
ax.hlines(1e-8, 0, 1)
plt.legend()
plt.yscale("log")
ax.set_xlim(0, 100)
ax.set_ylim(1e-9, 1e-5)
plt.grid(True)
ax.set_xlabel(r"Recovery efficiency (%)", fontsize=14)
ax.set_ylabel(r"[Complex] (M)", fontsize=14)
ax.set_title(
    r"pbSEC kinetic scheme, [P] = 10 $\mathrm{\mu}$M, [L] = 8 $\mathrm{\mu}$M, K$_\mathrm{D}$ = "
    + str(KD)
    + " µM",
    fontsize=16,
)
plt.show()
