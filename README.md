# pbSEC_simulation
Code repository for simulation of plate-based size exclusion chromatography (pbSEC) runs. 

This repository accompanies a soon to be published protocol.

![pbSEC simulation](https://raw.githubusercontent.com/stevenshave/pbSEC_simulation/main/pbSEC-simulation.png "Simulation of pbSEC kinetic scheme")


## pbSEC_equations.py
Contains functions for simulation of pbSEC, specificically:
- one_to_one_binding
    - Generic equation for simulation of 1:1 binding.
- one_to_one_binding_mpf
    - Generic equation for simulation of 1:1 binding taking high precision mpf parameters.
- pbSEC
    -  Simulate a round of pbSEC.
- pbSEC_simulate_n_rounds
    - Simulate N rounds of pbSEC (kinetic mode).
- pbSEC_iterate_until_undetectable
    - Perform multiple kinetic rounds of pbSEC until compound falls bellow a defined detection limit.

## 01_run_single_pbSEC.py
Simulates a simple pbSEC round and outputs compound concentrations

## 02_run_n_kinetic_rounds.py
Simulate multiple iterations of the kinetic pbSEC scheme

## 03_run_kinetic_rounds_until_undetectable.py
Simulate running rounds of kinetic pbSEC until compound falls beneath a specified detection limit.

## 04_plot_pbSEC_pKD_vs_complex.py and 04_plot_pbSEC_pKD_vs_logcomplex.py
Plot pKD vs complex concentration (non-log and log concentration versions)

## 05_plot_pbSEC_recovery_vs_logcomplex.py
Plot compound recovery rate vs logcomplex

## 06_plot_pbSEC_pKD_vs_num_rounds_detectable.py
Plot pKD vs number of rounds the compound is detectable for.

## 07_plot_pbSEC_3D_pKD_vs_logcomplex_vs_recovery_efficiency.py
3D plot of pKD vs recovery vs complex