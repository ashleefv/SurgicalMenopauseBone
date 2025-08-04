%% Solve model
% Solve model of ODEs for cell types, bone variables and chemicals.
% t- time (days)

% dependent variable definitions
% 1. PB - preosteoblasts
% 2. PC - preosteoclasts
% 3. C  - osteoclast
% 4. B  - osteoblast 
% 5. S  - osteocytes
% 6. Sc  - sclerositin
% 7. BMC  - bone mineralistion constant
% 8. Bd - bone density

% e- estrogen (imposed concentration as a function of time)

% solves ode for given initial condition, parameters and cases
% returns time, solution array, normalised solution array and estrogen vector

function [T,sol,sol_norm,BMD_norm,estrogen_vec]=solve_model(params,initialcond, t_array,t_ref,if_surgical,if_new_effects)

% ODE solver dynamic system with initial condition from fsolve.
options=odeset('RelTol',1e-13,'AbsTol',1e-13); % solver tolerances
[T,sol] = ode45(@(t,s) ode_rhs(t,s,params,if_surgical,if_new_effects),t_array,initialcond,options );

% calculate BMD at t=t_ref for normalisation 
BMD_ref=sol(T == t_ref,7)*params.BMC_0;
% normalise the solution vector relative to the initial condition
sol_norm=sol./initialcond;

% calculate the bone mineralisation density as product of bone min content and bone density
BMD_norm=sol(:,7).*params.BMC_0./BMD_ref;  % normalise by BMD at t_ref
estrogen_vec=zeros(size(T)); % empty storage vector
for i=1:length(T)
    % filll each array entry with estrogen value
    estrogen_vec(i)=estrogen(T(i),if_surgical,params.t_e,params.tau_e,params.k_dec,params.k_syn);
end
end


