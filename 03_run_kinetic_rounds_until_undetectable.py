"""
Simulate iterative pbSEC runs until instrument detection limit is reached

Perform iterative pbSEC runs using protein and ligand with known KD and iterate
until the specified detection limit is reached.
"""
from pbSEC_equations import pbSEC_iterate_until_undetectable

l = 8
kd = 43.5
p = 10
recovery_efficiency = 1.0
num_iterations = 10
l_detection_limit = 0.001

print(
    f"Using the following parameters:\n{p=}\n{l=}\n{kd=}\n{l_detection_limit=}\n{recovery_efficiency=},"
)
print(
    "Amount complex (also amount ligand) at each iteration =",
    pbSEC_iterate_until_undetectable(
        p, l, kd, l_detection_limit, recovery_efficiency=recovery_efficiency
    ),
)
