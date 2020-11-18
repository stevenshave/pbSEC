"""
Simulate a single pbSEC run

Perform a singular pbSEC run using protein and ligand with known KD
"""
from pbSEC_equations import pbSEC, one_to_one_binding

p = 10
l = 8
kd = 10
recovery_efficiency = 1.0


print(f"Using the following parameters:\n{p=}\n{l=}\n{kd=}\n{recovery_efficiency=}")
print(
    "Amount complex (also amount ligand) =",
    pbSEC(p, l, kd, recovery_efficiency=recovery_efficiency),
)
print(one_to_one_binding(10, 8, 10))
