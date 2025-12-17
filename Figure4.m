
% Generates Figure 4 for paper. 

clear all
close all

read_in_parameters
Spine_data%: load in menopause data.
if_new_effects=0; % turn off new effects so we can find steady states. 

set(groot,'DefaultLineMarkerSize',20)
% load in new params 
refit_bone_params=load('refit_bone_params.mat');%from 'Refit_bone_params.m'
refit_bone_params=refit_bone_params.refit_bone_params; % overwrite struct 

figure3=figure('units','inch','position',[0,0,7 ,6]);

BMD_ref=0.7606; % BMD at menopause onset. 
alpha=1; % set sclerostin levels to normal
params.e_PC=refit_bone_params(1);%
params.e_Sc=refit_bone_params(2);%
estrogen_level=1; % premenopause
initialcond = get_initial_condition_stab(params, estrogen_level,alpha);
Sc0=initialcond(6);
C0=initialcond(3);
B0=initialcond(4);

dBMD0=B0*params.lambda_B*repression(Sc0,params.sc_Omega)*(1+params.nu_Omega*activation(C0,params.r_Omega))-C0*params.lambda_C;
Ss=initialcond(5)*repression(estrogen_level,params.e_Sc)*alpha';

plot(100,100*dBMD0*365/BMD_ref,'k*',DisplayName='Premenopause',MarkerSize=20); hold on 

Ss0=Ss; % store this as refernce value for later plots


% now post menopause levels.
estrogen_level=params.k_syn/params.kappa_E;
initialcond = get_initial_condition_stab(params, estrogen_level,alpha);
Sc0=initialcond(6);
C0=initialcond(3);
B0=initialcond(4);
Ss=initialcond(5)*repression(estrogen_level,params.e_Sc)*alpha;
dBMD0=B0*params.lambda_B*repression(Sc0,params.sc_Omega)*(1+params.nu_Omega*activation(C0,params.r_Omega))-C0*params.lambda_C;
plot(Ss/Ss0*100,100*dBMD0*365/BMD_ref,'k.',DisplayName='Post-surgical menopause - no new effects',MarkerSize=30); hold on 



% Now create curve sweeping through alpha: varying sclerostin production
% rates
rates=[];
scs=[];
alphas=linspace(0.7,1.2,100);
for i=1:length(alphas)
estrogen_level=params.k_syn/params.kappa_E;

alpha=alphas(i);
initialcond = get_initial_condition_stab(params, estrogen_level,alpha);
Sc0=initialcond(6);
C0=initialcond(3);
B0=initialcond(4);
dBMD0=B0*params.lambda_B*repression(Sc0,params.sc_Omega)*(1+params.nu_Omega*activation(C0,params.r_Omega))-C0*params.lambda_C;
rates(i)=dBMD0*365;

Ss(i)=initialcond(5)*repression(estrogen_level,params.e_Sc)*alpha;
end

plot(Ss/Ss0*100,100*rates/BMD_ref,'k-',HandleVisibility='off'); hold on 

% Now plot 80% production 
alpha=0.8;
params.e_PC=refit_bone_params(1);%
params.e_Sc=refit_bone_params(2);%
estrogen_level=1;
initialcond = get_initial_condition_stab(params, estrogen_level,alpha);
Sc0=initialcond(6);
C0=initialcond(3);
B0=initialcond(4);

dBMD0=B0*params.lambda_B*repression(Sc0,params.sc_Omega)*(1+params.nu_Omega*activation(C0,params.r_Omega))-C0*params.lambda_C;
Ss=initialcond(5)*repression(estrogen_level,params.e_Sc)*alpha;

plot(Ss/Ss0*100,100*dBMD0*365/BMD_ref,'b.',DisplayName='Minimum production level - long term fit',MarkerSize=30); hold on 

% Now plot 94% production
alpha=0.94;
params.e_PC=refit_bone_params(1);%
params.e_Sc=refit_bone_params(2);%
estrogen_level=1;
initialcond = get_initial_condition_stab(params, estrogen_level,alpha);

Sc0=initialcond(6);
C0=initialcond(3);
B0=initialcond(4);

dBMD0=B0*params.lambda_B*repression(Sc0,params.sc_Omega)*(1+params.nu_Omega*activation(C0,params.r_Omega))-C0*params.lambda_C;
Ss=initialcond(5)*repression(estrogen_level,params.e_Sc)*alpha;

plot(Ss/Ss0*100,100*dBMD0*365/BMD_ref,'g.',DisplayName='Minmum production level - short term fit',MarkerSize=30); hold on 

lk=xline(80);
lk.LineWidth=2;
lk.LineStyle='-.';
lk.Color='b';
lk.HandleVisibility='off';

lk=xline(94);
lk.LineWidth=2;
lk.LineStyle='-.';
lk.Color='g';
lk.HandleVisibility='off';


legend(Location='southeast')
lk=yline(0);
lk.LineWidth=2;
lk.LineStyle='-.';
lk.Color=[0.5,0.5,0.5];
lk.DisplayName='No bone loss';

xlabel('Sclerostin production level (\% of premenopause level)')
ylabel('Steady BMD change (\% per year)')
exportgraphics(figure3,'Fig4.pdf', resolution = 300)


function initialcond = get_initial_condition_stab(params,estrogen_level,alpha)
% steady state solver without Bone mineral density
s0_fsolve = ones(1,6); % create vector to pass to fsolve
% optio
% solve for the steady state value of model without bone density when e=1.
options = optimset('Display','off');
x = fsolve(@(s) steady_state_rhs_stab(s,params,estrogen_level,alpha),s0_fsolve,options);

initialcond=[x,1]; %initial condition for ode
end

function s_prime= steady_state_rhs_stab(var,params,estrogen_level,alpha)
% RHS of model, removing bone density, can be solved for steady stateinitial condition,

% extract variables
var_cell=num2cell(var);
[PB,PC,C,B,S,Sc]=var_cell{:};

% Regulatory factors are proportional to osteoclast
r = C;
% initial constant estrogen conc.
E=estrogen_level;

t=0; % set time to be zero for time-dependent RHS
% RHS of odes


dPB = 1 - PB*repression(Sc,params.sc_PB)*params.omega_PB;
dPC = 1 - PC*repression(E,params.e_PC)*activation(Sc,params.sc_PC)*params.omega_PC;
dB = PB*repression(Sc,params.sc_PB)*params.omega_PB - (params.eta_B + params.omega_B)*B;
dSc = alpha*repression(E,params.e_Sc)*S- params.kappa_Sc*Sc; 
dC =  PC*params.omega_PC*repression(E,params.e_PC)*activation(Sc,params.sc_PC)- C*(params.eta_C);
dS = params.omega_B*B - params.eta_S*S; 

%combine to array
s_prime = [dPB dPC dC dB dS dSc]';

end