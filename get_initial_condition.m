%s_prime = [dpb dpc dc db ds dsc,dBd]';
%[pb,pc,c,b,s,sc,BMC,Bd,scb]=var_cell{:};
% 1. PB - preosteoblasts
% 2. PC - preosteoclasts
% 3. C  - osteoclasts
% 4. B  - osteoblasts 
% 5. S  - osteocytes
% 6. Sc  - sclerostin
% 7. Bd - bone density
function initialcond = get_initial_condition(params, if_new_effects)
% steady state solver without Bone mineral density
s0_fsolve = ones(1,6); % create vector to pass to fsolve
% optio
% solve for the steady state value of model without bone density when e=1.
options = optimset('Display','off');
x = fsolve(@(s) steady_state_rhs(s,params,if_new_effects),s0_fsolve,options);

initialcond=[x,1]; %initial condition for ode
end

function s_prime= steady_state_rhs(var,params,if_new_effects)
% RHS of model, removing bone density, can be solved for steady stateinitial condition,

% extract variables
var_cell=num2cell(var);
[PB,PC,C,B,S,Sc]=var_cell{:};

% Regulatory factors are proportional to osteoclast
r = C;
% initial constant estrogen conc.
e = 1;


t=0; % set time to be zero for time-dependent RHS
% RHS of odes
RHS_of_equations

%combine to array
s_prime = [dPB dPC dC dB dS dSc]';

end






