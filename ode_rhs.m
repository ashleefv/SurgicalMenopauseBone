% RHS of ode for time dependent solver
function s_prime= ode_rhs(t,var,params,if_surgical,if_new_effects)

% extract variables
var_cell=num2cell(var);
[PB,PC,C,B,S,Sc,Bd]=var_cell{:};

% Regulatory factors are proportional to osteoclast
r = C;

% time dependent estrogen conc.
e=estrogen(t,if_surgical,params.t_m,params.tau_e,params.k_deg,params.k_syn);


% RHS of odes
RHS_of_equations 

%combine to array
s_prime = [dPB dPC dC dB dS dSc,dBd]';

end

