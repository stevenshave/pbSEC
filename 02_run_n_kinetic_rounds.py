"""
Simulate iterative pbSEC runs

Perform iterative pbSEC runs using protein and ligand with known KD
"""
from pbsec import pbSEC_simulate_n_rounds

l = 8
kd = 15.2
p = 10
recovery_efficiency = 1.0
num_iterations = 10

print(
    f"Using the following parameters:\n{p=}\n{l=}\n{kd=}\n{num_iterations=}\n{recovery_efficiency=},"
)
print(
    "Amount complex (also amount ligand) at each iteration =",
    pbSEC_simulate_n_rounds(
        p, l, kd, num_iterations, recovery_efficiency=recovery_efficiency
    ),
)
