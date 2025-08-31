% Human Lumbar spine BMD data
% % Novemeber 2023 10, 2023


% This script loads in time and BMD values for all experimental papers used
% for fitting to

% Initialization
% read in parameters and other settings

read_in_parameters;
t_surgery = params.t_m/365; % surgery time in years
t_ref=params.t_m; % Normalisation BMD value.

% Data import

% ALL TIME DATA IN days ALL bmd DATA NORMALISED.

% Pansini et al 1995
% Human data

% Pansini monitored for 12 years following since menses or surgery
t_pansini=[0,1,2,3,4,12]*365; %+ t_surgery*365;  %time in days since surgical meno

% timepoints are taken as the end point of binned time categories in the pansini paper 
 
% Percentages of BMD at years
% surgical menopause
BMD_pansini=[1,1-0.059,1-0.134,1-0.185,1-.215,1-.258];
BMD_pansini_SD=[0,0.117,0.111,0.096,0.082,0.134];

% Natural Menopause
BMD_N_pansini=[1,1-.104,1-.146,1-.153,1-.168,1-.218]; % natural menopause 
BMD_N_pansini_nat_SD=[0,0.134,0.142,0.166,0.127,0.147]; 

% -----------------------------------------------------------------
 % Hilber et al. 2016 surgical menopause human
% % Ranges are available in paper

t_hibler=[0,18]*30; % time in days since surgical menopause
BMD_hibler=[1,1-0.09]; %  Normalised on baseline amount at surgery

% SD calculated from CI with N=52, CI=[1-0.12,1-0.07] SD=0.03*sqrt(52)/1.95
BMD_hibler_SD=[0,0.03*sqrt(52)/1.95];% SD not listed for baseline

%------------------------------------------------------------
% Hajidakis 2002 (measured the women at the same age and also listed their onset
% time, didnt track the women over time)

t_hajidakis=[47.5,52.5,57.5,62.5,67.5]*365; % middle of age categories measured at

%mean onset ages for each category (in years)
t_hajidakis_onset_N=49.1;
t_hajidakis_onset_E=38.3;
t_hajidakis_onset_S=38.1;

% Normalise all data on the amount of BMD in the natural case at 47.5
% BMD_norm=8.022222222222222;

BMD_norm=8.022222222222222*0.1; %(normalise on amount as in Yasui 2007 perimenopausal.)

% Natural meno
BMD_N_hajidakis=[8.022222222222222, 7.533333333333331,6.955555555555554, 6.666666666666666, 6.2]*0.1./BMD_norm;
mean_SD=[9.777777777777775, 8.666666666666664, 7.999999999999998,...
    8.044444444444443, 7.288888888888888]*0.1./BMD_norm;
BMD_N_hajidakis_SD=(mean_SD-BMD_N_hajidakis);;

% Surgical

BMD_hajidakis=[6.755555555555554, 6.755555555555554, 6.711111111111111, 7.244444444444442, 7.155555555555555
]*0.1./BMD_norm;
mean_SD=[ 7.5555555555555545, 7.7333333333333325, 7.577777777777777, 8.622222222222222, 8.777777777777777
]*0.1./BMD_norm;
BMD_hajidakis_SD=(mean_SD-BMD_hajidakis);


%early onset
BMD_E_hajidakis=[6.666666666666666, 6.333333333333332, 6.444444444444444,...
    6.422222222222222, 6.444444444444444]*0.1./BMD_norm;
mean_SD=[ 7.7333333333333325,7.155555555555555, 7.5555555555555545...
,7.399999999999998, 7.866666666666665]*0.1./BMD_norm;
BMD_E_hajidakis_SD=(mean_SD-BMD_hajidakis);

%---------------------------------------------------------------------------
% Ohta 2002 
BMD_norm=1.059;


% natural 
t_N_ohta=[50.99 58.79]*365; % age of women when measured
t_since_menses_N=[13.9,100]*30; % time since menses for natural in two measurement cases (Days)
t_ohta_onset_N=t_N_ohta-t_since_menses_N;

% surgical 
t_ohta=[47.99 55.99]*365; % age of women when measured
t_since_menses=[12.9,112]*30;  % time since surgery in two measurement cases (Days)
t_ohta_onset=t_ohta-t_since_menses;


% natural 
BMD_N_ohta=[0.961,0.826]./BMD_norm;
BMD_N_ohta_SD=[0.123,0.122]./BMD_norm;

% surgical 
BMD_ohta=[0.931,0.7978102189781022]./BMD_norm;
BMD_ohta_SD=[0.115, 0.0723]./BMD_norm;
 % final point extracted using image analysis 

 
% Interval after menopause/
% oophorectomy (months)
% 
% 13.99/3.6 100.99/39.4 12.99/8.6 112.29/89.5

%---------------------------------------------------------------------------
% Looker 1998: From Jorg, hip BMD used to calibrate Jorg model
% modifed from 'looker1998_and_recknor2015_placebo.csv' (additional columns removed, rows with no BMD value removed)
% Time (days),BMD total hip,Error low: BMD total hip,Error high: BMD total hip,P1NP,Error low: P1NP,Error high: P1NP,CTX,Error low: CTX,Error high: CTX,Medication: Blosozumab (mg)

% Think they took this data from Looker 1998 Table 7 (Standardized total
% femur)
data_looker_import=[9125,1.0
12775,0.9895287958115183
16425,0.9633507853403143
20075,0.9172774869109946
23725,0.8471204188481676
23893,0.8374267877301257
24089,0.8420591424252776
24173,0.8381988468417487
24257,0.8376841407616858
24341,0.8376841407616858
24453,0.8358826694984084];


data_looker=data_looker_import(1:5,:); % keep only the first 5 entries which are pre-treatment, from Fig. 2 (green dots)
SD_looker = [123, 130, 136, 139, 140]; 

% Only take the last two points (after menopause -- which is not included)
t_looker = [data_looker(1:5,1)]'; 
BMD_norm = 955;  
t_looker_onset = 50; %taken from eyeballing Jorg

 
BMD_looker=data_looker(1:5,2)/0.959297; % normalised by BMD at 50 in jorg
mean_SD_looker=SD_looker(1:5)./BMD_norm;
BMD_looker_SD=mean_SD_looker';

%---------------------------------------------------------------------------
% Yasui 2007: Longitudinal study with 21 bilaterally oophrectomized women
% (measured BMD from L2-L4 vertebrate, PTH, estrone and estradiol

% Measurements (post-op): 1 month, 6 months, 12 months, 2-5 years, 6-10
% years
t_yasui_surgery=[48.3, 48.3, 47.5, 48.1, 50.3 ]*365; % middle of age categories measured at
t_yasui_pre_surg= [48.1, 48.0, 46.5, 45.4, 41.6]*365; 
% Normalise all data on the amount of BMD premenopausal
BMD_norm=1.059;
t_yasui_pre =  47.9; 
BMD_yasui=[0.999, 1.003, 0.967, 0.909, 0.8]./BMD_norm;
mean_SD=[ 0.087, 0.092, 0.129, 0.139, 0.087]./BMD_norm;
BMD_yasui_SD=(mean_SD-BMD_yasui);

%---------------------------------------------------------------------------
% Chittacharoen 1997: 50 surg, 50 perimenopausal, no hormonal treatment,
%

% Normalise all data on the amount of BMD premenopausal
BMD_norm=1.153;
t_chitt_pre =  48.98*365; 
t_chitt_surg = 50.16*365; 
% Measurements (post-op): <= 3 years,   4-6 years, 7-9 years, 10-12, >=12 
% years (taken from Table 3)    
t_chitt_surgery_delta=[0, 3, 5, 8, 11, 13]*365; % middle of age categories measured at
BMD_chitt=[BMD_norm, 1.155, 0.991, 1.016, 0.888, 0.893]./BMD_norm;
mean_SD_chitt=[0.149, 0.094, 0.096, 0.134, 0.163, 0.123]./BMD_norm;



%---------------------------------------------------------------------------
% Chittacharoen 1999: 309 natural 102 surgical, free of metabolic disease,
% duration of menopause on average is 5.18 \pm 5.12, 5.67 \pm 4.76

% Only gives range in a table? Not time points...


%---------------------------------------------------------------------------
% Yasui 2007: Longitudinal study with 21 bilaterally oophrectomized women
% (measured BMD from L2-L4 vertebrate, PTH, estrone and estradiol

% Measurements (post-op): 1 month, 6 months, 12 months, 2-5 years, 6-10
% years
t_yasui_surgery=[48.3, 48.3, 47.5, 48.1, 50.3 ]*365; % middle of age categories measured at
t_yasui_pre_surg= [48.1, 48.0, 46.5, 45.4, 41.6]*365; 
% Normalise all data on the amount of BMD premenopausal
BMD_norm=1.059;
t_yasui_pre =  47.9; 
BMD_yasui=[0.999, 1.003, 0.967, 0.909, 0.8]./BMD_norm;
mean_SD=[ 0.087, 0.092, 0.129, 0.139, 0.087]./BMD_norm;
BMD_yasui_SD=(mean_SD-BMD_yasui);

%---------------------------------------------------------------------------
% Chittacharoen 1997: 50 surg, 50 perimenopausal, no hormonal treatment,
%

% Normalise all data on the amount of BMD premenopausal
BMD_norm=1.153;
t_chitt_pre =  48.98*365; 
t_chitt_surg = 50.16*365; 
% Measurements (post-op): <= 3 years,   4-6 years, 7-9 years, 10-12, >=12 
% years (taken from Table 3)    
t_chitt_surgery_delta=[0, 3, 5, 8, 11, 13]*365; % middle of age categories measured at
BMD_chitt=[BMD_norm, 1.155, 0.991, 1.016, 0.888, 0.893]./BMD_norm;
mean_SD_chitt=[0.149, 0.094, 0.096, 0.134, 0.163, 0.123]./BMD_norm;

