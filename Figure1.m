
clear all 
close all
% Generates Figure 1 for paper. 

Spine_data; % Import data 

% Combine into 1 array: Natural data 
t_N=[t_pansini/365,t_looker(4:5)/365-50 t_hajidakis/365-t_hajidakis_onset_N,t_N_ohta/365-t_N_ohta_onset/365];
BMD_N=[BMD_N_pansini, BMD_looker(4:5)', BMD_N_hajidakis,BMD_N_ohta]*100;

% Combine into 1 array: Surgical 
t_S=[t_pansini/365, t_hibler/365,t_hajidakis/365-t_hajidakis_onset_S,...
    t_yasui_surgery/365-t_yasui_pre_surg/365, t_chitt_surgery_delta/365,...
    t_ohta/365-t_ohta_onset/365];
BMD_S=[BMD_pansini, BMD_hibler,BMD_hajidakis ,BMD_yasui,BMD_chitt, BMD_ohta]*100;


% Calulcate linear fit for first 15 years of BMD dynamics for each dataset.
% fit a line through 100 both datasets 

t_N_exd=t_N(t_N<=15);
BMD_N_exd=BMD_N(t_N<=15);
fun = @(a) sum((BMD_N_exd-(a*t_N_exd+100)).^2);% Natural dataset fitting
[xN, ~, ~, ~] = fminunc(fun,-0.1);
"the percentage decrease in natural meno in BMD in 15 years is ", xN

t_S_exd=t_S(t_S<=15);
BMD_S_exd=BMD_S(t_S<=15);

funs_exd = @(a) sum((BMD_S_exd-(a*t_S_exd+100)).^2);
 
[xS_exd, fval, exitflag, output] = fminunc(funs_exd,-0.1);
"the percentage decrease in surgical meno in BMD in 15 years is ", xS_exd

% Plot colours and cmap for subplots (b,c)
cmap = hot(13);

pansini_color       = cmap(1,:); 
pansini_marker      = 'o';
hajidakis_color     = cmap(2,:); 
hajidakis_marker    = '^'; 
ohta_color          = cmap(3,:); 
ohta_marker         = '^'; 
looker_color        = cmap(4,:); 
looker_marker       = 'o'; 
hibler_color        = cmap(5,:); 
hibler_marker       = 'square'; 
yasui_color         = cmap(6,:); 
yasui_marker        = 'v'; 
chitt_color         = cmap(7,:); 
chitt_marker        = 'diamond'; 

set(groot,'DefaultAxesFontSize',19);




% Plot data (Natual dotted, surgical solid, early dash-dots)
close all
figure_1=figure('units','inch','position',[0,0,16,6]);
tt = tiledlayout(1,3); 




% t=tiledlayout(1,2);
ax2=nexttile(tt);

% First tile: Normalised on time since last menses data used for estimation

errorbar(t_pansini/365, BMD_N_pansini*100 ,BMD_N_pansini_nat_SD*100,...
       'marker',pansini_marker,'linestyle','none','MarkerSize',10,'MarkerFaceColor',pansini_color, 'color',pansini_color, 'linewidth',2,'DisplayName','Pansini et al.'); hold on


errorbar(t_hajidakis/365-t_hajidakis_onset_N, BMD_N_hajidakis*100 ,BMD_N_hajidakis_SD*100,...
   'marker',hajidakis_marker,'markerfacecolor',hajidakis_color,'linestyle','none','color',hajidakis_color,'MarkerSize',10, 'linewidth',2,'DisplayName','Hajidakis et al.');
t=linspace(0,15,100);

plot(t,xN*t+100,'k-.','DisplayName',['Fit: ', num2str(round(xN,2)), '\%/yr']); hold on

errorbar(t_N_ohta/365-t_N_ohta_onset/365, BMD_N_ohta*100,BMD_N_ohta_SD*100,...
    'marker',ohta_marker, 'linestyle','none','color', ohta_color,...
   'MarkerSize',10, 'linewidth',2,'DisplayName','Ohta et al.');


errorbar(t_looker/365-t_looker_onset, BMD_looker*100,BMD_looker_SD*100,...
   'marker',looker_marker,'linestyle','none', 'color',looker_color,...
   'MarkerSize',10, 'linewidth',2,'DisplayName','Looker et al.');

ylim([60,140])
xlim([-5,30])


title('(a)')
xlabel('Years since menopause onset')
ylabel('Normalised BMD \%')

leg=legend;
title(leg,'Natural menopause')
leg.Location='northeast';
leg.NumColumns = 2;



ax3=nexttile(tt);

% First tile: Normalised on time since last menses data used for estimation
errorbar(t_pansini/365, BMD_pansini*100,BMD_pansini_SD*100,'marker',pansini_marker,'linestyle','none','MarkerFaceColor',pansini_color, 'color',pansini_color,'MarkerSize',10,'linewidth',2,...
    'DisplayName', 'Pansini et al.'); hold on

errorbar(t_hibler/365, BMD_hibler*100 ,BMD_hibler_SD*100,...
   'marker',hibler_marker,'linestyle','none','color',hibler_color,'MarkerSize',10, 'linewidth',2,'DisplayName','Hibler et al.');


errorbar(t_ohta/365-t_ohta_onset/365, BMD_ohta*100, BMD_ohta_SD*100,...
   'marker',ohta_marker,'linestyle','none','color',ohta_color,'MarkerSize',10, 'linewidth',2,'DisplayName','Ohta et al.');

errorbar(t_yasui_surgery/365-t_yasui_pre_surg/365, BMD_yasui*100,mean_SD*100,...
   'marker',yasui_marker,'markerfacecolor',yasui_color,'linestyle','none','color',yasui_color,	'MarkerSize',10, 'linewidth',2,'DisplayName','Yasui et al.');


errorbar(t_hajidakis/365-t_hajidakis_onset_S, BMD_hajidakis*100 ,BMD_hajidakis_SD*100,...
    'marker',hajidakis_marker,'markerfacecolor',hajidakis_color,'linestyle','none','color',hajidakis_color,'MarkerSize',10, 'linewidth',2,'DisplayName','Hajidakis et al.');


 errorbar(t_chitt_surgery_delta/365, BMD_chitt*100,mean_SD_chitt*100,...
   'marker',chitt_marker,'markerfacecolor','none','linestyle','none','color',chitt_color,	'MarkerSize',10, 'linewidth',2,'DisplayName','Chittacharoen et al.');


plot(t,xS_exd*t+100,'k-.','DisplayName',['Fit: ', num2str(round(xS_exd,2)), '\%/yr'])
title('(b)')
xlabel('Years since menopause onset')
ylabel('Normalised BMD \%')



ylim([60,140])
xlim([-5,30])


leg=legend;
title(leg,'Surgical menopause')
leg.Location='north';
leg.NumColumns = 2; 
% Create the primary axis in the tile
axPrimary = nexttile(tt);
box(axPrimary, 'on');
axis(axPrimary, 'equal'); 
% Generate estrogen concentrations in both natural and surgical cases. 
t=linspace(0,57*365,1000); % 57 years
% t=linspace(params.t_m-10,params.t_m+10,1000);
es=zeros(size(t));
en=zeros(size(t));
for i=1:length(t)
    es(i)=estrogen(t(i),1,params.t_m,params.tau_E,params.kappa_E,params.k_syn);
        en(i)=estrogen(t(i),0,params.t_m,params.tau_E,params.kappa_E,params.k_syn);
end
%
% es=estrogen(t,1,params.t_m,params.tau_E,params.kappa_E,params.k_syn);
% % en=estrogen(t,0,params.t_m,params.tau_E,params.kappa_E,params.k_syn);

plot(axPrimary, t/365-params.t_m/365,en,'r','LineWidth',2); hold on
plot(axPrimary,t/365-params.t_m/365,es,'k','LineWidth',2)
title({'(c)'})
xlim([-2,30])
% Choose section to isolate
xSection = [0, 25];  
ySection = [20, 35];  

xlabel('Years since menopause onset')
ylabel('Relative blood estrogen conc.')

leg=legend('Natural menopause','Surgical menopause');
leg.Location="north";
axSecondary = axes('Position',[0.8 0.5 0.13 0.22], 'Box', 'on'); 
axis(axSecondary, 'equal');

% Zoom in to the 10 days around onset

t=linspace(params.t_m-10,params.t_m+10,1000);
es=zeros(size(t));
en=zeros(size(t));
for i=1:length(t)
    es(i)=estrogen(t(i),1,params.t_m,params.tau_E,params.kappa_E,params.k_syn);
        en(i)=estrogen(t(i),0,params.t_m,params.tau_E,params.kappa_E,params.k_syn);
end
%
plot(axSecondary, t-params.t_m,en,'r',t-params.t_m,es,'k','LineWidth',2)
xlim([-1,3])
xlabel('Days since menopause onset')
title('')





tt.TileSpacing = 'compact';
tt.Padding = 'compact';

set(gcf,'color','w')
fig = figure(1); 
fig.Position = [0.4167 0.8611 17.6389 6.8194]; 




% save figure manually and move legend panel 3 right
exportgraphics(figure_1,'Fig1.pdf',Resolution = 300)

