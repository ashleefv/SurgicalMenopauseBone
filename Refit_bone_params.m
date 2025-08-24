% Script estimating parameters in Jorg model that are related to BMD
% dynamics. Use 4 data sets of natural menopause (including Looker)

% Target parameters: lambda_B - formation rate of bone, lambda_c - bone resorption rate,
% e_PC - Preosteoclast suceptibility to estrogen, e_Sc, Sclerostin
% susceptibiity to estrogen

clear all
close all; 
read_in_parameters
tend = 80*365; %final time
tstart=0*365; %initial time
t_ref=21*365;  % calibration time 25 years

Spine_data%: load in menopause data.

% Combine into 1 array: Natural
t_N=[t_pansini/365, t_looker/365-50,t_hajidakis/365-t_hajidakis_onset_N,t_ohta/365-t_ohta_onset/365];
BMD_N=[BMD_N_pansini, BMD_looker', BMD_N_hajidakis,BMD_ohta]*100;

% Combine into 1 array: Surgical 
t_S=[t_pansini/365, t_hibler/365,t_hajidakis/365-t_hajidakis_onset_S,...
    t_yasui_surgery/365-t_yasui_pre_surg/365, t_chitt_surgery_delta/365,...
    t_N_ohta/365-t_ohta_onset_N/365];
BMD_S=[BMD_pansini, BMD_hibler,BMD_hajidakis ,BMD_yasui,BMD_chitt, BMD_N_ohta]*100;

%%------------ Panel a data: nat and surg no new effects  --------
params.omega_ovx=0;
params.tau=0;
params.eta_ovx=0;
params.t_m=27*365;
if_new_effects=0;
if_surgical=0;
initialcond = get_initial_condition(params,if_new_effects);

%%% Estimate parameters:  parameters related to estrogen action and
%%% sclerostin

lb=[0, 0,0,0];  %lower bound of parameters  


% Initial guess to nonlinear least squares
% lambda_B     = 1.29E-06; % ND, Threshold for estrogen action on osteoclasts
% lambda_C    = 3.82E-06; % ND, Threshold for estrogen action on pre-osteoclasts
% [lambda_B, lambda_C, e_PC, e_Sc]
kguess=[1e-6, 4e-7, 0.5, 1]; 

%implementing lsqcurvefit: look for optimal params

[sorted_t_N_vector,sorting] = sort(t_N*365); 
sorted_BMD_N_vector = BMD_N(sorting); 


OPTIONS = optimoptions('lsqnonlin','StepTolerance',1e-12,...
 'FunctionTolerance',1e-16,'optimalitytolerance', 1e-12,...
 'MaxFunctionEvaluations', 200,'Algorithm','levenberg-marquardt',...
 'Display','iter');


% Sorted t data vector is in YEARS not DAYS
[refit_bone_params,resnorm,residual,exitflag,output] = lsqnonlin(@(k) ...
     solve_model_varyparam_k(k,params,initialcond, tstart:1:tend,t_ref, sorted_t_N_vector, if_surgical,...
   if_new_effects, sorted_BMD_N_vector/100),...
   kguess, lb, [Inf,Inf,Inf,Inf], OPTIONS)

% save data for use in later scripts.
save('refit_bone_params.mat',"refit_bone_params",'-mat')

% Function which ode for given search parameters, returns error to data

function [F]=solve_model_varyparam_k(k,params,initialcond, ...
    t_array,t_ref,t_data,if_surgical,if_new_effects,  ...
     BMD_N_vector)...


params.lambda_B = k(1); 
params.lambda_C = k(2); 
params.e_PC = k(3); 
params.e_Sc = k(4); 


% ODE solver dynamic system with initial condition from fsolve.
options=odeset('RelTol',1e-13,'AbsTol',1e-13); % solver tolerances
[T,sol] = ode45(@(t,s) ode_rhs(t,s,params,if_surgical,if_new_effects),t_array,initialcond,options );

% calculate BMD at t=25 for normalisation as in Jorg
BMD_25=sol(T == t_ref,7)*params.BMC_0; 

%Calculate solution at the data points

[T,sol] = ode45(@(t,s) ode_rhs(t,s,params,if_surgical,if_new_effects),t_data,initialcond,options );

% normalise the solution vector relative to the initial condition
sol_norm=sol./initialcond;

% calculate the bone mineralisation density as product of bone min content and bone density
BMD_norm=sol(:,7).*params.BMC_0./BMD_25;  % normalise by BMD at age 25.
BMD_norm = BMD_norm'; 

F = BMD_norm - BMD_N_vector; 
% create vector of estrogen value at time points in solution
estrogen_vec=zeros(size(T)); % empty storage vector
for i=1:length(T)
    % filll each array entry with estrogen value
    estrogen_vec(i)=estrogen(T(i),if_surgical,params.t_m,params.tau_E,params.kappa_E,params.k_syn);
end


end





