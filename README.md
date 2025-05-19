# SurgicalMenopauseBone

Analysis scripts for paper "Mathematical Modeling of Bone Remodeling in Surgical Menopause"


# Authors
Anna C. Nelson(1), Edwina Yeo(2)
(1) Department of Mathematics \& Statistics, University of New Mexico, Albuquerque, NM, USA
(2)Department of Mathematics, University College London, United Kingdom


# Overview
To run all analysis and create all figures evaluate main.m. This repositiory contains an ODE
based model which simulates the dynamics of bone cells in the BMU and predits BMD density over time.

# Data import files
a. Spine_data.m: contains surgical and natural menopause data from humans studies of spinal BMD
b. read_in_parameters.m reads model parameters2 and sets plotting parameters
c. params.txt: read in by read_in_parameters.m

# Functions
 a. get_initial_condition.m: solves model at steady state and returns normalisation values for each cells species
  b. solve_model.m: solves model returning cell species and BMD percentage over time. 
 c. activation.m & repression.m are hill functions 
 d. estrogen.m returns either relative natural or surgical estrogen concs over time.
 e. ode_rhs.m is the ODE function passed to solve_model.m
 f. RHS_of_equations.m generates the RHS of the equations to pass to ode_rhs and get_initial_condition

# 3. Analysis scripts
   a. Refit-bone-params.m: re-estimates bone parameters using wide natural menopause dataset
  b. Fit_new_params.m: estimates parameters for new effects using surgical menopause datasets.
 c. main.m: runs all analysis and creates figures. 

# 4. Plotting scripts:
  a. Figure1.m and Figure2.m generate paper figures. 
