% Script estimating parameters in Jorg model that are related to model
% response to estrogen. 
% dynamics. Use 4 data sets of natural menopause (including Looker)

% Target parameters: % e_PC - Preosteoclast suceptibility to estrogen, e_Sc, Sclerostin
% susceptibiity to estrogen

clear all
% close all; 
read_in_parameters
Spine_data%: load in menopause data.

tend = 50*365; %final time
tstart=-30*365; %initial time
t_ref=0*365-0.001;  % calibration time 0 years- shifted by a small amount to avoid coincinding with datapoints
params.t_m=0*365;% menopause  time 0 years

% Combine into 1 array: Natural
t_N=[t_pansini/365, t_looker/365-50,t_hajidakis/365-t_hajidakis_onset_N,t_N_ohta/365-t_N_ohta_onset/365];
BMD_N=[BMD_N_pansini, BMD_looker', BMD_N_hajidakis,BMD_N_ohta]*100;

% Combine into 1 array: Surgical 
t_S=[t_pansini/365, t_hibler/365,t_hajidakis/365-t_hajidakis_onset_S,...
    t_yasui_surgery/365-t_yasui_pre_surg/365, t_chitt_surgery_delta/365,...
    t_ohta/365-t_ohta_onset/365];
BMD_S=[BMD_pansini, BMD_hibler,BMD_hajidakis ,BMD_yasui,BMD_chitt, BMD_ohta]*100;
[sorted_t_N_vector,sorting] = sort(t_N*365); 
sorted_BMD_N_vector = BMD_N(sorting); 

%%------------ Panel a data: nat and surg no new effects  --------
params.omega_ovx=0;
params.tau=0;
params.eta_ovx=0;
if_new_effects=0;
if_surgical=0;

initialcond = get_initial_condition(params,if_new_effects);
%%% Estimate parameters:  parameters related to estrogen action and
% Initial guess to nonlinear least squares
% [ e_PC, e_Sc]
kguess=[0.9377, 9.5954];  % start at the values of Jorg et al.

%implementing lsqcurvefit: look for optimal params
OPTIONS = optimoptions('lsqnonlin','StepTolerance',1e-8,...
 'FunctionTolerance',1e-8,'optimalitytolerance', 1e-8,...
 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000,'Algorithm','levenberg-marquardt');

t_solve=sort([tstart,t_ref,sorted_t_N_vector,tend]); 
[refit_bone_params,resnorm,residual,exitflag,output] = lsqnonlin(@(k) ...
     solve_model_varyparam_k(k,params,initialcond, t_solve,t_ref, sorted_t_N_vector/365, if_surgical,...
   if_new_effects,sorted_BMD_N_vector/100),...
   kguess,[0.05,0.05],[20,20],OPTIONS);

% save data for use in later scripts.
save('refit_bone_params.mat',"refit_bone_params",'-mat')




function [F]=solve_model_varyparam_k(k,params,initialcond, ...
    t_array,t_ref,t_data,if_surgical,if_new_effects,  ...
     BMD_N_vector)...


params.e_PC = k(1); 
params.e_Sc = k(2); 


initialcond = get_initial_condition(params,if_new_effects); % solve for initial condition of model

% ODE solver dynamic system with initial condition from fsolve.
options=odeset('RelTol',1e-8,'AbsTol',1e-8); % solver tolerances
[T,sol] = ode45(@(t,s) ode_rhs(t,s,params,if_surgical,if_new_effects),t_array,initialcond,options );

% calculate BMD at t=25 for normalisation as in Jorg
BMD_25=sol(T == t_ref,7)*params.BMC_0;
BMD_norm=sol(:,7).*params.BMC_0./BMD_25(1);  % normalise by BMD at age 25.
BMD_norm = BMD_norm';

% create a mask so we only compare the model at the data points
mask =  t_array == t_ref;
mask(1)=1;
mask(end)=1;

F = (BMD_norm(mask<1)  - BMD_N_vector);



end





