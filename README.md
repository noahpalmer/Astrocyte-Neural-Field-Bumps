This repository contains the code needed to generate the results for the paper "Astrocytic resource diffusion stabilizes persistent activity in neural fields" by Noah Palmer, Heather Cihak, Daniele Avitabile and Zachary Kilpatrick (to be submitted).

**Overview:**
The code within this repository can be used to simulate and analyze stationary bump solutions in an astrocyte-neural field model. It includes scripts for finite difference simulation of the astrocyte-neural field model, perturbation analysis, and numerical approximation of stability boundaries in phase space.

**Files:**
- astrocyte_neural_field_sim.m :
  Contains the finite difference scheme used to verify the existence of stationary bump solutions as well as run simulations with a small perturbation.
- numericalphase.m :
  Contains the code necessary to numerically simulate small perturbations to stationary bump solutions in various   parameter regimes. Returns bump velocity and drift.
- analyticphase.m :
  Contains code to numerically estimate the boundary between unstable and stbale regions for both finite        diffusion and the large diffusion limit, using the Evans function derived in Sections 4.5 and 4.6.1 of the paper.
- fourierphase.m :
  Estimates the stability boundaries for finite diffusion using a low-Fourier mode truncation method outlined in Section 5.3.
- detcomp.m :
  Auxiliary function used to compute the Evans function for analyticphase.m.

**Requirements:**
All code was run on MATLAB Version: 24.1.0.2603908 (R2024a) Update 3.

**Simulation length:**
The simulations carried out in "numericalphase.m" can take a long time depending on resolution. On a laptop running a Windows 11 Home operating system a 50 by 50 resolution simulation takes approximately 
