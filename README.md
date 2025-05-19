# SurgicalMenopauseBone

Analysis scripts for paper "Mathematical Modeling of Bone Remodeling in Surgical Menopause"


# Authors
Anna C. Nelson(1), Edwina Yeo(2), Ashlee N. Ford Versypt (3)
- (1) Department of Mathematics \& Statistics, University of New Mexico, Albuquerque, NM, USA
- (2) Department of Mathematics, University College London, United Kingdom
- (3) Department of Biomedical Engineering, University at Buffalo, The State University of New York, Buffalo, NY, USA


# Overview
To run all analysis and create all figures evaluate main.m. This repositiory contains an ODE
based model which simulates the dynamics of bone cells in the BMU and predits BMD density over time.

# Data import files
- Spine_data.m: contains surgical and natural menopause data from humans studies of spinal BMD
- read_in_parameters.m reads model parameters2 and sets plotting parameters
- params.txt: read in by read_in_parameters.m

# Functions
 - get_initial_condition.m: solves model at steady state and returns normalisation values for each cells species
 - solve_model.m: solves model returning cell species and BMD percentage over time.
 - activation.m & repression.m are hill functions
 - estrogen.m returns either relative natural or surgical estrogen concs over time.
 - ode_rhs.m is the ODE function passed to solve_model.m
 - RHS_of_equations.m generates the RHS of the equations to pass to ode_rhs and get_initial_condition

# Analysis scripts
- Refit-bone-params.m: re-estimates bone parameters using wide natural menopause dataset
- Fit_new_params.m: estimates parameters for new effects using surgical menopause datasets.
- main.m: runs all analysis and creates figures. 

# Plotting scripts:
- Figure1.m and Figure2.m generate paper figures. 
