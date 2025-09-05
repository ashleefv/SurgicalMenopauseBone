%% Read in parameters
% Read in parameters from "Params.txt" 
% Save all parameters in a struct "params"

% Open text file and extract string names and parameter values
fileID1 = fopen('params.txt','r');
c1 = textscan(fileID1,'%s = %f; %*[^\n]');
C = [c1{1},num2cell(c1{2})].';
params = struct(C{:});


params.t_m =0*365;      % days, Rescaled time parameter, onset time in days: 0

% Estrogen decay surgical menopause parameters
E_0 = 156;  % picogram/ml initial concentration of estrogen
kappa_E = log(2)*(24*60)/161; % days^(-1), Rescaled timescale of estrogen decline, estogen half life
% 0.065 is a rounded version of log(2)*15/161. Slighlty more precise version of k_syn is 15/e_0*log(2)*(24*60)/161
k_syn = 0.065*(24*60)/E_0; % days^(-1), Rescaled timescale 
params.E_0 = E_0;

% terms for estrogen effect on apoptosis
Eovx=k_syn/kappa_E; % post ovx estrogen level;
params.Eovx = Eovx;
params.kappa_E = kappa_E;
params.k_syn = k_syn;

%parameters for new effects set to zero.
params.eta_ovx=0;  % maximum amount of increase in apoptosis after surgery
params.tau=0; % rate at which cyte apoptosis returns to normal: days^(-1)
params.omega_ovx=0;  % maximum amount of increase in apoptosis after surgery


% Plotting parameters: sets all plots to have these properties
set(groot,'DefaultTextInterpreter','latex') 
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultAxesFontSize',16);
set(groot,'DefaultLineLineWidth',2)