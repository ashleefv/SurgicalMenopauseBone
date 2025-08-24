%% Right hand side of the euqations
% This file defines the RHS of the surgical menopause model

% PB pre-osteoblasts
% PC pre-osteoclasts
% C osteoclasts
% B osteoblasts
% Sc sclerostin
% BMC bone mineral content
% Bd Bone density


dPB = 1 - PB*repression(Sc,params.sc_PB)*params.omega_PB;
dPC = 1 - PC*repression(E,params.e_PC)*activation(Sc,params.sc_PC)*params.omega_PC;
dB = PB*repression(Sc,params.sc_PB)*params.omega_PB - (params.eta_B + params.omega_B)*B;
dSc = repression(E,params.e_Sc)*S- params.kappa_Sc*Sc; 
dBd=B*params.lambda_B*repression(Sc,params.sc_Omega)*(1+params.nu_Omega*activation(r,params.r_Omega))-C*params.lambda_C;
dC =  PC*params.omega_PC*repression(E,params.e_PC)*activation(Sc,params.sc_PC)- C*(params.eta_C);
dS = params.omega_B*B - params.eta_S*S; 


% If new effects are on then use different functional forms 
if if_new_effects

    % apoptosis rate:
    eta_surg=params.eta_S.*(t<params.t_m)...
        +((params.eta_S+params.eta_S*params.eta_ovx...
        *exp(-params.tau*abs(t-params.t_m)))).*(t>=params.t_m);

    % overwrite osteocyte equation
    dS = params.omega_B*B - eta_surg.*S;

    % Increased differentation of PC 
     params.omega_surg=params.omega_PC.*(t<params.t_m)...
        +((params.omega_PC+params.omega_ovx...
        *exp(-params.tau*abs(t-params.t_m)))).*(t>=params.t_m);

    % overwrite preosteoclasts and osteoclast equation
    dPC   = 1 - PC*repression(E,params.e_PC)*activation(Sc,params.sc_PC)*params.omega_surg;
    dC =  PC*params.omega_surg*repression(E,params.e_PC)*activation(Sc,params.sc_PC)- C*(params.eta_C);
        
end



