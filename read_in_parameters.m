%% Read in parameters
% Read in parameters from "Params.txt" 
% Save all parameters in a struct "params"

% Open text file and extract string names and parameter values
fileID1 = fopen('params.txt','r');
c1 = textscan(fileID1,'%s = %f; %*[^\n]');
C = [c1{1},num2cell(c1{2})].';
params = struct(C{:});


params.t_e =27*365;      % days, Rescaled time parameter, onset time in days: 38 years

% estrogen decay surgical menopause parameters
e_0 = 156;  % picogram/ml initial concentration of estrogen
k_dec = log(2)*(24*60)/161; % days^(-1), Rescaled timescale of estrogen decline, estogen half life
k_syn = 0.065*(24*60)/e_0; % days^(-1), Rescaled timescale 
params.e_0 = e_0;

% terms for estrogen effect on apoptosis
eovx=k_syn/k_dec; % post ovx estrogen level;
params.eovx = eovx;
params.k_dec = k_dec;
params.k_syn = k_syn;

%parameters for new effects set to zero.
params.eta_ovx=0;  % maximum amount of increase in apoptosis after surgery
params.tau=0; % rate at which cyte apoptosis returns to normal: days^(-1)
params.omega_ovx=0;  % maximum amount of increase in apoptosis after surgery

% external constant rate of sclerostin production: needs to be set to the value of S
% outside of steady state:
params.omega_Sc=1; 

% Plotting parameters: sets all plots to have these properties
set(groot,'DefaultTextInterpreter','latex') 
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultAxesFontSize',16);
set(groot,'DefaultLineLineWidth',2)