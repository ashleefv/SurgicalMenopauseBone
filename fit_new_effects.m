% Script estimating parameters new model changes using surgical menopause data

% Target parameters: 
% params.eta_ovx  % maximum amount of increase in osteocyte apoptosis after surgery
% % params.tau=; % rate at which osteocyte apoptosis & increased differentiation returns to normal: days^(-1)
% params.omega_ovx  % maximum amount of increase in preosteoclast differentation after surgery

clear all
close all; 

read_in_parameters % load model parameters
tend = 80*365; %final time for simulation
tstart=0*365; %initial time
t_ref=21*365;  % calibration time 21 years

Spine_data % Load in menpause data

% Combine into 1 array, surgical data
t_S=[t_pansini/365, t_hibler/365,t_hajidakis/365-t_hajidakis_onset_S,...
    t_yasui_surgery/365-t_yasui_pre_surg/365, t_chitt_surgery_delta/365,...
    t_N_ohta/365-t_ohta_onset_N/365];
BMD_S=[BMD_pansini, BMD_hibler,BMD_hajidakis ,BMD_yasui,BMD_chitt, BMD_N_ohta]*100;

% Use fitted bone parameters from 'Refit-bone-params.m'

refit_bone_params=load('refit_bone_params.mat');%from 'Refit_bone_params.m'
refit_bone_params=refit_bone_params.refit_bone_params; % overwrite struct 

% Case 2: natural menpause with new parameters 

params.lambda_B=refit_bone_params(1); %
params.lambda_C=refit_bone_params(2);%
params.e_PC=refit_bone_params(3);%
params.e_Sc=refit_bone_params(4);%
if_surgical=1; % turn on new estrogen dynamics 
if_new_effects=1; % turn on new effects

% Case 1: fitting to all data 
kguess=[10,0.001,1]; % Intial guess for all params
lb=[0,1e-5,0]; % lower bound
ub=[1e2,1,1e2];% upper bbestound
[sorted_t_S_vector,sorting] = sort(t_S*365); 
sorted_BMD_S_vector = BMD_S(sorting); 
tsolve=[linspace(tstart,t_ref-1,10),t_ref,linspace(t_ref+1,tend,10)];
initialcond = get_initial_condition(params,if_new_effects);

 OPTIONS = optimoptions('lsqcurvefit','StepTolerance',1e-16,...
 'FunctionTolerance',1e-14,'optimalitytolerance', 1e-14,'MaxFunctionEvaluations',4000,'Algorithm','levenberg-marquardt');
 % OPTIONS = optimoptions('lsqcurvefit','StepTolerance',1e-8,...
 % 'FunctionTolerance',1e-8,'optimalitytolerance', 1e-8,'MaxIterations',5000,'MaxFunctionEvaluations',2000,'Algorithm','levenberg-marquardt');


[best_params_all,fit_error_long,~,~] = lsqcurvefit(@(k,tdata) ...
    new_params_err(k,params,initialcond, tsolve,t_ref,  sorted_t_S_vector/365, if_surgical,...
   if_new_effects) , kguess,   sorted_t_S_vector/365,sorted_BMD_S_vector/100,lb,ub,OPTIONS);
"FIT error long is "
fit_error_long


% Case 2: fitting to data up to 15 years 
sorted_BMD_S_vector=sorted_BMD_S_vector(sorted_t_S_vector/365<15);
sorted_t_S_vector=sorted_t_S_vector(sorted_t_S_vector/365<15);

initialcond = get_initial_condition(params,if_new_effects);
[best_params_short,fit_error_short,~,~] = lsqcurvefit(@(k,tdata) ...
    new_params_err(k,params,initialcond, tstart:1:tend,t_ref,  sorted_t_S_vector/365, if_surgical,...
   if_new_effects) , kguess,   sorted_t_S_vector/365,sorted_BMD_S_vector/100,lb,ub,OPTIONS);

"FIT error short is "
fit_error_short

% save data for use in later scripts.
save('fit_effects_long.mat',"best_params_all",'-mat')
save('fit_effects_short.mat',"best_params_short",'-mat')

function [F]=new_params_err(k,params,initialcond, t_array,t_ref,t_data,if_surgical,if_new_effects)

params.eta_ovx = k(1); 
params.tau = k(2); 
params.omega_ovx=k(3);

% ODE solver dynamic system with initial condition from fsolve.
options=odeset('RelTol',1e-13,'AbsTol',1e-13); % solver tolerances
[T,sol] = ode45(@(t,s) ode_rhs(t,s,params,if_surgical,if_new_effects),t_array,initialcond,options );
% 
% calculate BMD at t=21 for normalisation
BMD_25=sol(T == t_ref,7)*sol(T == t_ref,8);
BMD_norm=sol(:,7).*sol(:,8)./BMD_25;  % normalise by BMD at age 21.
BMD_norm = BMD_norm';

% BMD_25=0.6185;
% plot(T/365-params.t_e/365,BMD_norm,'r-.','DisplayName','new fit'); hold on
% plot(T,BMD_norm,'r-.','DisplayName','new fit'); hold on

%Calculate solution at the data points
t_data_scaled=params.t_e+t_data*365;



[T,sol] = ode45(@(t,s) ode_rhs(t,s,params,if_surgical,if_new_effects),[0,t_ref,t_data_scaled],initialcond,options );

% normalise the solution vector relative to the initial condition
sol_norm=sol./initialcond;

% calculate the bone mineralisation density as product of bone min content and bone density
BMD_norm=sol(:,7).*sol(:,8)./BMD_25;  % normalise by BMD at age 25.
BMD_norm = BMD_norm';


F = BMD_norm(3:end);% - BMD_S_vector/100

% create vector of estrogen value at time points in solution
estrogen_vec=zeros(size(T)); % empty storage vector
for i=1:length(T)
    % filll each array entry with estrogen value
    estrogen_vec(i)=estrogen(T(i),if_surgical,params.t_e,params.tau_e,params.k_dec,params.k_syn);
end

end



