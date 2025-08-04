%% Right hand side of the euqations
% This file defines the RHS of the surgical menopause model

% PB pre-osteoblasts
% PC pre-osteoclasts
% C osteoclasts
% B osteoblasts
% Sc sclerostin
% BMC bone mineral content
% Bd Bone density


dPB = 1 - PB*repression(Sc,params.s_PB)*params.omega_PB;
dPC = 1 - PC*repression(e,params.e_PC)*activation(Sc,params.s_PC)*params.omega_PC;
dB = PB*repression(Sc,params.s_PB)*params.omega_PB - (params.eta_B + params.omega_B)*B;
dSc = repression(e,params.e_Sc)*S- params.kappa_S*Sc; 
dBd=B*params.lambda_B*repression(Sc,params.s_Omega)*(1+params.nu_Omega*activation(r,params.r_Omega))-C*params.lambda_C;
dC =  PC*params.omega_PC*repression(e,params.e_PC)*activation(Sc,params.s_PC)- C*(params.eta_C);
dS = params.omega_B*B - params.eta_S*S; 


% If new effects are on then use different functional forms 
if if_new_effects

    % apoptosis rate:
    eta_surg=params.eta_S.*(t<params.t_e)...
        +((params.eta_S+params.eta_S*params.eta_ovx...
        *exp(-params.tau*abs(t-params.t_e)))).*(t>=params.t_e);

    % overwrite osteocyte equation
    dS = params.omega_B*B - eta_surg.*S;

    % Increased differentation of PC 
     params.omega_surg=params.omega_PC.*(t<params.t_e)...
        +((params.omega_PC+params.omega_ovx...
        *exp(-params.tau*abs(t-params.t_e)))).*(t>=params.t_e);

    % overwrite preosteoclasts and osteoclast equation
    dPC   = 1 - PC*repression(e,params.e_PC)*activation(Sc,params.s_PC)*params.omega_surg;
    dC =  PC*params.omega_surg*repression(e,params.e_PC)*activation(Sc,params.s_PC)- C*(params.eta_C);
        
end



