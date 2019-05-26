# Water_and_Energy_budgets
Enforce the closure of the water and energy budgets simultaneously at every grid

This script enforces the closure of the surface water and energy budgets simultaneously.
Given a set of imbalanced components of the water and energy budgets (net radiation Rn, sensible heat flux H, latent heat flux/evapotranspiration LH, ground heat flux G, precipitation P, runoff Q and change in water storage DS) and their associated uncertainties, the function Solve_Budgets_Grid adjusts each of these component based on their relative uncertainties, then returns a suite of optimized (adjusted) fluxes and optimized uncertainties that solve the two budgets simultaneously.
This script was used to derive a new dataset Conserving Land-Atmosphere Synthesis Suite (CLASS). The algorithm is explained in:
Hobeichi, S., Abramowitz, G., Evans, J. P.: Conserving Land-Atmosphere Synthesis Suite (CLASS): Journal of Climate, In review, 2019.
