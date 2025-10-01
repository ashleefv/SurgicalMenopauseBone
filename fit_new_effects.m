% Script estimating parameters new model changes using surgical menopause data

% Target parameters: 
% params.eta_ovx  % maximum amount of increase in osteocyte apoptosis after surgery
% params.tau=; % rate at which osteocyte apoptosis & increased differentiation returns to normal: days^(-1)
% params.omega_ovx  % maximum amount of increase in preosteoclast differentation after surgery

clear all
close all; 
Spine_data % Load in menpause data

read_in_parameters % load model parameters

tend = 50*365; %final time
tstart=-30*365; %initial time
t_ref=0*365-1e-3;  % calibration time 0 years
params.t_m=0*365;% calibration time 0 years


% Combine into 1 array, surgical data
t_S=[t_pansini/365, t_hibler/365,t_hajidakis/365-t_hajidakis_onset_S,...
    t_yasui_surgery/365-t_yasui_pre_surg/365, t_chitt_surgery_delta/365,...
    t_ohta/365-t_ohta_onset/365];
BMD_S=[BMD_pansini, BMD_hibler,BMD_hajidakis ,BMD_yasui,BMD_chitt, BMD_ohta]*100;

% Use fitted bone parameters from 'Refit-bone-params.m'

refit_bone_params=load('refit_bone_params.mat');%from 'Refit_bone_params.m'
refit_bone_params=refit_bone_params.refit_bone_params; % overwrite struct 

% Case 2: natural menpause with new parameters 

params.e_PC=refit_bone_params(1);%
params.e_Sc=refit_bone_params(2);%
if_surgical=1; % turn on new estrogen dynamics 
if_new_effects=1; % turn on new effects


% Case 1: fitting to all data 
kguess=[0.5,0.001,0.5]; % Intial guess for all params
lb=[0,0,0]; % lower bound: [no effect, permanent effect, no effect]
ub=[5,1,5];% upper bound: [5 times increased apoptosis, effect which lasts 1 day, 5 times increase apop]

[sorted_t_S_vector,sorting] = sort(t_S*365); 
sorted_BMD_S_vector = BMD_S(sorting); 

t_solve_S=sort([tstart,t_ref,sorted_t_S_vector,tend]);


OPTIONS = optimoptions('lsqnonlin','StepTolerance',1e-8,...
 'FunctionTolerance',1e-8,'optimalitytolerance', 1e-8,...
 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000,'Algorithm','levenberg-marquardt');

[best_params_all,resnorm,residual,exitflag,output]= lsqnonlin(@(k) ...
    new_params_err(k,params, t_solve_S,t_ref,  if_surgical,...
   if_new_effects,sorted_BMD_S_vector/100) , kguess, lb,ub,OPTIONS)
RMSE_long = sqrt(resnorm/length(sorted_BMD_S_vector))*100 % units of BMD percentage
%
save('fit_effects_long.mat',"best_params_all",'-mat')
%
% % Case 2: fitting to data up to 15 years 
sorted_BMD_S_vector=sorted_BMD_S_vector(sorted_t_S_vector/365<15);
sorted_t_S_vector=sorted_t_S_vector(sorted_t_S_vector/365<15);
t_solve_S=sort([tstart,t_ref,sorted_t_S_vector,tend]);

% initialcond = get_initial_condition(params,if_new_effects);
[best_params_short,fit_error_short,~,~] = lsqnonlin(@(k) ...
    new_params_err(k,params, t_solve_S,t_ref,  if_surgical,...
   if_new_effects,sorted_BMD_S_vector/100) , kguess, lb,ub,OPTIONS);

RMSE_short = sqrt(fit_error_short/length(sorted_BMD_S_vector))*100 % units of BMD percentage

save('fit_effects_short.mat',"best_params_short",'-mat')

function [F]=new_params_err(k,params, t_array,t_ref,if_surgical,if_new_effects, BMD_vector)

params.eta_ovx = k(1); 
params.tau = k(2); 
params.omega_ovx=k(3);

initialcond = get_initial_condition(params,if_new_effects); % solve for initial condition of model

% ODE solver dynamic system with initial condition from fsolve.
options=odeset('RelTol',1e-8,'AbsTol',1e-8); % solver tolerances

[T,sol] = ode45(@(t,s) ode_rhs(t,s,params,if_surgical,if_new_effects),t_array,initialcond,options );

% calculate BMD at t=25 for normalisation as in Jorg
BMD_25=sol(T == t_ref,7)*params.BMC_0;
BMD_norm=sol(:,7).*params.BMC_0./BMD_25;  % normalise by BMD at age 25.
BMD_norm = BMD_norm';

% o we only compare the model at the data points
mask =  t_array == t_ref;
mask(1)=1;
mask(end)=1;
F = (BMD_norm(mask<1)  - BMD_vector);
end



