% Master script which runs all analysis files and creates figures. 

tic
% 1. Re-estimate bone params
Refit_bone_params
toc

tic
% 2. Estimate new model params
%fit_new_effects
toc
%
% 3. Generate figs
Figure1
Figure3
