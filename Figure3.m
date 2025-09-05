
% Generates Figure 3 for paper. 

clear all
close all; 
read_in_parameters; % load model parameters
Spine_data % load menopause data


tend = 50*365; %final time
tstart=-30*365; %initial time
t_ref=0*365;  % calibration time 0 years
params.t_m=0*365;% menopause time 0 years



% Combine into one array for plotting: Natural data without the Looker et
% al dataset.
t_N=[t_pansini/365, t_hajidakis/365-t_hajidakis_onset_N,t_ohta/365-t_ohta_onset/365]; % years since onset
BMD_N=[BMD_N_pansini, BMD_N_hajidakis,BMD_ohta]*100; % Relative BMD percentage


% Combine into one array for plotting: Surgical data 
t_S=[t_pansini/365, t_hibler/365,t_hajidakis/365-t_hajidakis_onset_S,...
    t_yasui_surgery/365-t_yasui_pre_surg/365, t_chitt_surgery_delta/365,...
    t_N_ohta/365-t_ohta_onset_N/365];% years since onset
BMD_S=[BMD_pansini, BMD_hibler,BMD_hajidakis ,BMD_yasui,BMD_chitt, BMD_N_ohta]*100;% Relative BMD percentage


% Generate model solutions ------------------------------------
% ---------------------------------------------------------------
%

%Data for plot 2a-----------
% Case 1: natural menopause with Jorg et al. parameters.
% solve models and plot
params.omega_ovx=0; % turn off new effects - no increased differentiation
if_new_effects=0; % turn off new effects - no increased apoptosis
if_surgical=0; % natural meno
initialcond = get_initial_condition(params,if_new_effects); % solve for initial condition of model
% solve for BMD dynamics
[T_n_pre,~,~,BMD_norm_pre,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);


% save data for use in later scripts.
refit_bone_params=load('refit_bone_params.mat');%from 'Refit_bone_params.m'
refit_bone_params=refit_bone_params.refit_bone_params; % overwrite struct 
% save data for use in later scripts.
best_params_all=load('fit_effects_long.mat'); %from 'Fit_new_effects.m'
best_params_short=load('fit_effects_short.mat'); %from 'Fit_new_effects.m'
best_params_all=best_params_all.best_params_all ;% overwrite struct 
best_params_short=best_params_short.best_params_short ;% overwrite struct 

% Case 2: natural menpause with new parameters 

% params.lambda_B=refit_bone_params(1); %
% params.lambda_C=refit_bone_params(2);%
params.e_PC=refit_bone_params(1);%
params.e_Sc=refit_bone_params(2);%

initialcond = get_initial_condition(params,if_new_effects);
[T_n,~,sol_norm,BMD_normn,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);

% case 3: surgical menopause using new bone parameters but no new effects. 

if_surgical=1; % turns on surgical menopause
initialcond = get_initial_condition(params,if_new_effects);
[T_s,~,sol_norm_fit,BMD_norms,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);
%

%Data for plot 2b-----------

% Long fit surgical data with new effects
if_new_effects=1; 
if_surgical=1;

params.eta_ovx =best_params_all(1); % Increased percentage of apoptosis
params.tau = best_params_all(2); % timescale of effects
params.omega_ovx=best_params_all(3);
initialcond = get_initial_condition(params,if_new_effects);

[~,~,sol_fit_long,BMD_fit_long,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);

% Create upper bound by increasing or decreasing effects by 25%: increase
% percentages and increase duration by decreasing timescale by 25%.
params.eta_ovx =best_params_all(1)*1.25; 
params.tau = best_params_all(2)*0.75; 
params.omega_ovx=best_params_all(3)*1.25;
initialcond = get_initial_condition(params,if_new_effects);
[~,~,~,BMD_fit_long_upper,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);

% Create lower bound by increasing or decreasing effects by 25%: decrease
% percentages and decrease duration by increasing timescale by 25%.
params.eta_ovx =best_params_all(1)*0.75; 
params.tau = best_params_all(2)*1.25; 
params.omega_ovx=best_params_all(3)*0.75;
initialcond = get_initial_condition(params,if_new_effects);

[~,~,~,BMD_fit_long_lower,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);


% Short fit surgical data with new effects
params.eta_ovx =best_params_short(1); 
params.tau = best_params_short(2); 
params.omega_ovx=best_params_short(3);
initialcond = get_initial_condition(params,if_new_effects);

[~,~,sol_fit_short,BMD_fit_short,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);

 % Create upper bound by increasing or decreasing effects by 25%: increase
% percentages and increase duration by decreasing timescale by 25%.
params.eta_ovx =best_params_short(1)*1.25; 
params.tau = best_params_short(2)*0.75; 
params.omega_ovx=best_params_short(3)*1.25;
initialcond = get_initial_condition(params,if_new_effects);

[~,~,~,BMD_short_upper,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);

% Create lower bound by increasing or decreasing effects by 25%: decrease
% percentages and decrease duration by increasing timescale by 25%.
params.eta_ovx =best_params_short(1)*0.75; 
params.tau = best_params_short(2)*1.25; 
params.omega_ovx=best_params_short(3)*0.75;
initialcond = get_initial_condition(params,if_new_effects);

[~,~,~,BMD_short_lower,~]=solve_model(params, initialcond, tstart:1:tend, t_ref, if_surgical, if_new_effects);

%
% Create plot-----------------------------------------------------------

close all

figure3=figure('units','inch','position',[0,0,16,6]);
t1 = tiledlayout(1,3,'TileSpacing','Compact');
t2 = tiledlayout(t1,'flow','TileSpacing','Compact');
t3 = tiledlayout(t1,'flow','TileSpacing','Compact');
t4 = tiledlayout(t1,'flow','TileSpacing','Compact');
t4.Layout.Tile = 3;
t3.Layout.Tile = 2;
t2.Layout.Tile = 1;


% Plot a
nexttile(t2);
plot(t_S,BMD_S,'k.','MarkerSize',25,'DisplayName','SM data'); hold on
plot(t_looker/365-50, BMD_looker'*100,'ro','DisplayName','NM data Looker et al.',MarkerSize=10)
plot(t_N,BMD_N,'r.','MarkerSize',25,'DisplayName','NM data other sources'); hold on
plot(T_n/365-params.t_m/365,BMD_norm_pre*100,'r-.','DisplayName','Model NM:  Jorg et al.')
plot(T_n/365-params.t_m/365,BMD_normn*100,'r','DisplayName','Model NM: refit Jorg et al.')
plot(T_s/365-params.t_m/365,BMD_norms*100,'k','DisplayName','Model SM: no new effects')
ylim([60,120])
xlim([-5,30])
legend()

xlabel('Years since menopause onset')
ylabel('Relative BMD \%')
title("(a)")


%Plot b
nexttile(t3)

long_BMD=BMD_S(t_S>15);
long_t=t_S(t_S>15);



plot(t_S,BMD_S,'k.','MarkerSize',25,'DisplayName','SM data'); hold on
plot(long_t,long_BMD,'b.','MarkerSize',25,'DisplayName','Long term SM data'); hold on
plot(T_s/365-params.t_m/365,BMD_norms*100,'k','DisplayName','Model SM: no new effects')
plot(T_s/365-params.t_m/365,BMD_fit_short*100,'g','DisplayName','Model SM:  fit up to 15 years')

plot(T_s/365-params.t_m/365,BMD_fit_long*100,'b','DisplayName','Model SM: fit up to 30 years')

xlim([-5,30])
ylim([60,120])

kk=patch([T_s(1:200:end);flip(T_s(1:200:end));T_s(1)]/365-params.t_m/365,...
    [BMD_fit_long_lower(1:200:end);flip(BMD_fit_long_upper(1:200:end));BMD_fit_long_lower(1)]*100,'b',FaceAlpha=0.1,EdgeColor='none');
kk.HandleVisibility='off';

kk=patch([T_s(1:200:end);flip(T_s(1:200:end));T_s(1)]/365-params.t_m/365,...
    [BMD_short_upper(1:200:end);flip(BMD_short_lower(1:200:end));BMD_short_upper(1)]*100,'g',FaceAlpha=0.1,EdgeColor='none');
kk.HandleVisibility='off';
legend()

xlabel('Years since menopause onset')
ylabel('Relative BMD \%')
title("(b)")


%Plot ci
nexttile(t4)
plot(T_s/365-params.t_m/365,sol_norm_fit(:,3)*100,'k','DisplayName','Model SM: no new effects'); hold on
plot(T_s/365-params.t_m/365,sol_fit_short(:,3)*100,'g','DisplayName','Model SM:  fit up to 15 years'); hold on
plot(T_s/365-params.t_m/365,sol_fit_long(:,3)*100,'b','DisplayName','Model SM: fit up to 30 years')
legend()
ylabel('Osteoclasts \%')
title("(c)")
xlim([-5,30])

%Plot cii
nexttile(t4)
plot(T_s/365-params.t_m/365,sol_norm_fit(:,4)*100,'k'); hold on
plot(T_s/365-params.t_m/365,sol_fit_short(:,4)*100,'g')
plot(T_s/365-params.t_m/365,sol_fit_long(:,4)*100,'b'); hold on


ylabel('Osteoblast \%')
xlim([-5,30])

%Plot ciii
nexttile(t4)
plot(T_s/365-params.t_m/365,sol_norm_fit(:,5)*100,'k'); hold on

plot(T_s/365-params.t_m/365,sol_fit_short(:,5)*100,'g')
plot(T_s/365-params.t_m/365,sol_fit_long(:,5)*100,'b'); hold on


ylabel('Osteocytes \%')
xlabel('Years since menopause onset')
xlim([-5,30])

exportgraphics(figure3,'Fig3.pdf')
