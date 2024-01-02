"""
Plot pbSEC rounds, pKd vs num rounds detectable for

"""

import numpy as np
from pbsec import pbSEC_iterate_until_undetectable
import matplotlib.pyplot as plt

XAXIS_BEGINNING = 3  # pKD of 3 is 1 mM
XAXIS_END = 6  # pKD of 12 is 1 pM
NUM_POINTS_ON_XAXIS = 2000  # Publication used 2000 pts along X
PROTEIN_CONC = 10.0
LIGAND_CONC = 8.0  # Singular compound conc in pool
RECOVERY_EFFICIENCIES = np.arange(1.0, 0.4, -0.1)
L_DETECTION_LIMIT = 0.01


x_axis = np.linspace(XAXIS_BEGINNING, XAXIS_END, NUM_POINTS_ON_XAXIS)
ligand_kd_range = 10 ** (-x_axis) * 1e6
num_iterations_taken = np.full(
    (len(RECOVERY_EFFICIENCIES), NUM_POINTS_ON_XAXIS), np.nan, dtype=np.int
)

for recovery_efficiency_i, recovery_efficiency in enumerate(RECOVERY_EFFICIENCIES):
    for kd_i, kd in enumerate(ligand_kd_range):
        num_iterations_taken[recovery_efficiency_i, kd_i] = len(
            pbSEC_iterate_until_undetectable(
                PROTEIN_CONC,
                LIGAND_CONC,
                kd,
                L_DETECTION_LIMIT,
                recovery_efficiency=recovery_efficiency,
            )
        )

fig, ax = plt.subplots(1, 1, figsize=(7.204724, 5.09424929292))

for recovery_efficiency_i, recovery_efficiency in enumerate(RECOVERY_EFFICIENCIES):
    ax.plot(
        x_axis,
        num_iterations_taken[recovery_efficiency_i],
        label=f"{recovery_efficiency*100:2.0f} % recovery",
    )
plt.legend()

ax.set_xticklabels(
    ["3 (mM)", " ", "4", " ", "5", " ", r"6 ($\mathrm{\mu}$M)"]
)
ax.set_xlim(XAXIS_BEGINNING, XAXIS_END)
ax.set_ylim(0, 10)
plt.grid(True)
ax.set_xlabel(r"Ligand pK$_\mathrm{D}$", fontsize=16)
ax.set_ylabel(r"Num rounds detectable (N)", fontsize=16)
ax.set_title(
    "pbSEC kinetic scheme, [Prot] = 10 $\mathrm{\mu}$M, [Ligand] = 8 $\mathrm{\mu}$M,\n detection limit ="
    + str(L_DETECTION_LIMIT)
    + " $\mathrm{\mu}$M"
)
plt.show()
