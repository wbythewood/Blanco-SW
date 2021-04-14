% Script to setup parameters used for ASWMS

% prepend functions directory to MATLAB path
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
addpath(functionspath);

%IdString = 'TaTest_wide';
%IdString = 'AGU_7D';
IdString = 'MacTest';

% Set up paths
parameters.workingdir = '/Users/whawley/Research/Blanco-SW/';
parameters.figdir = [parameters.workingdir,'figures/ASWMS/',IdString,'/'];  % wbh addition - separate dir for figures
parameters.configDir = [parameters.workingdir,'config/']; %this is where the config files live
parameters.dataDir = [parameters.workingdir,'data/']; %data here
parameters.SacDbDir = [parameters.dataDir,'CORRSEIS_SAC/']; %the corrected sac files
%parameters.SacDbDir = [parameters.dataDir,'SAC_Events_TA/'];
%parameters.SacDbDir = [parameters.dataDir,'SAC_Events/'];
%parameters.MatDbDir = [parameters.dataDir,'eventmat/']; 
parameters.MatDbDir = [parameters.dataDir,'eventmat_',IdString,'/']; %the matlab version of the corrected sac files
parameters.ASWMSDir = [parameters.workingdir,'03_ASWMS/']; %everything after converting sac 2 mat goes in ASWMS dir
parameters.MatFilesDir = [parameters.ASWMSDir,'matfiles_',IdString,'/']; %modified mat files will go here
%parameters.MatFilesDir = [parameters.ASWMSDir,'matfiles_New/'];
parameters.MapsDir = '/Users/whawley/data/maps/'; %where some of my general map files are

% filenames
parameters.PACStaFile = [parameters.configDir,'stalist_PAC_nw.txt'];
parameters.JDFStaFile = [parameters.configDir,'stalist_JDF_nw.txt'];
parameters.StaFile = [parameters.configDir,'stalist_nw.txt'];
parameters.BadStaFile = [parameters.configDir,'stalist_bad.txt'];

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
%%%% Global settings
parameters.proj_name = 'Blanco';
parameters.component = 'LHZ';   % determined by filenames
parameters.lalim = [41.5 45.5];%[42.5 45] ;
parameters.lolim = [-131.5 -125.25];%[-131.5 -125.5];
%TA
%parameters.lalim = [33.0 38.0];%[42.5 45] ;
%parameters.lolim = [-120.0 -114.0];%[-131.5 -125.5];

parameters.gridsize = 0.25;   % in degrees
parameters.periods = [20 25 32 40 50 60 80 100 120 130 150]; 
%parameters.periods = [20 25 32 40 50 60 80 100];  % SP in seconds
%parameters.periods = [80 100 120 140 160]; %LP
% parameters.periods = round(logspace(log10(20),log10(150),15));
parameters.minSta = 5; % if fewer than this no of stations, evt skipped.

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% % parameters for data downloading (if using IRIS DMC)
% parameters.start_time = '2012-11-01 00:00:00';
% parameters.end_time = '2012-11-20 00:00:00'; % put '' for using 4 days before current date
% parameters.is_use_timestamp = 0;
% parameters.network = 'X9';
% parameters.minMw = 6.5;
% parameters.maxdepth = 500;
% parameters.datalength = 7200;  % in second
% parameters.resample_delta = 1; % in second
% parameters.dbpath = './sacdata/';
% parameters.eventfile = 'eventlist';

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% Parameters for own data selection criteria
parameters.dbpath = parameters.SacDbDir;
parameters.eventfile = [parameters.configDir,'BlancoEventList_M6.5.txt'];
%parameters.eventfile = [parameters.configDir,'BlancoProblematicEvt.txt'];
%parameters.eventfile = [parameters.configDir,'BlancoEventTest_LP.txt'];
parameters.minMw = 5.0;
parameters.maxdepth = 500;
parameters.snr_tol = 3;

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for the auto_win_select.m
parameters.dbpath = parameters.SacDbDir;
parameters.largest_epidist_range = 3000;
parameters.cycle_before = 2;
parameters.cycle_after = 3;
parameters.min_dist_tol = deg2km(20);
%parameters.max_dist_tol = deg2km(100); 
parameters.max_dist_tol = deg2km(160);
parameters.min_groupv = 2;
parameters.max_groupv = 5;
parameters.cent_freq = 0.025;
parameters.min_sta_num = 5; %10 %JBR

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for the cross-correlation measurement
% (gsdfmain.m)
parameters.minstadist = 5;
%parameters.maxstadist = 250; %200;   % station cross-correlation distance in km
parameters.maxstadist = 600;   % station cross-correlation distance in km
parameters.is_rm_resp = 0;
parameters.periods = sort(parameters.periods);  % make sure periods are ascending
parameters.refv = 4;   % to select the correct cycle
%parameters.refv = 5;   % to select the correct cycle
parameters.refphv = ones(size(parameters.periods))*4;
%parameters.refphv = ones(size(parameters.periods))*5;
parameters.min_width = 0.06;  % to build up gaussian filters
parameters.max_width = 0.10;  
parameters.wintaperlength = 30;   % taper to build up the isolation filter
parameters.prefilter = [10,160]; %[15,160]; %[10,200];
%parameters.prefilter = [10,110]; %SP
%parameters.prefilter = [70,200]; %LP
parameters.xcor_win_halflength = 300; %200;  % window for the cross-correlation
parameters.xcor_win_iter = zeros(size(parameters.periods)); % re-apply the xcor window due to measured group delay, should be same length as periods, not used anymore
parameters.Nfit = 2; %4; % 2
parameters.Ncircle = 5;
parameters.cohere_tol = 0.65; % minimum coherenecy between two stations
parameters.tp_tol = 10;  % seconds away from averaged phase velocity 
parameters.tp_tol = 20;  % seconds away from averaged phase velocity 

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for the tomography
% (eikonal_eq.m helmholtz_eq.m)
% parameters.smweight_array = 3*[0.4 0.3 0.2 0.2 0.2 0.5 1 2];  % smoothing weight for the deltaSx and delta Sy
parameters.smweight_array = 3*[0.4 0.3 0.2 0.2 0.2 0.5 1 2 2 3 3]; 
%parameters.smweight_array = 3*[0.4 0.3 0.2 0.2 0.2 0.5 1 2]; %SP
%parameters.smweight_array = 3*[1 2 2 3 3]; %LP
% parameters.smweight_array = 3*[0.4 0.3 0.2 0.2 0.2 0.5 0.5 0.5 1 1 1 2 3 3 3]; 
parameters.flweight_array = 100*ones(length(parameters.periods)); % JBR
parameters.flweight_array = 0.5*ones(length(parameters.periods)); % JBR
parameters.raydensetol=deg2km(parameters.gridsize)*2;
parameters.Tdumpweight = 0;  % dumping the ray to the girgle circle path
parameters.Rdumpweight = 0;  % dumping the region to have the same phase velocity
parameters.fiterrtol = 6; %3;   % error allowed in the wavelet fitting
parameters.isRsmooth = 1;  % smoothing due to Sx and Sy or Sr and S_theta
parameters.dterrtol = 4; %2;    % largest variance of the inversion error allowed
parameters.inverse_err_tol = 4; %2;  % count be number of standard devition
parameters.min_amp_tol = 0.1;  % station with amplitude smaller than this ratio of average amplitude will not be used.
parameters.amp_var_tol = 4; %2; % how much times change of amplitude of single station to the mean value of nearby stations should be considered as bad measurement
parameters.alpha_range = [1 1];
parameters.alpha_search_grid = 0.1;
parameters.min_stadist_wavelength = 0; %min wavelengths allowed between stations
parameters.max_stadist_wavelength = 999; %  max wavelengths
% azimuthal anisotropy details
parameters.aziweight = 0; %global weight
parameters.smweight0_azi = 1; % sm weight for azi
parameters.flweight0_azi = 0; % fl weight for azi

%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameter for stacking 
% (stack_phv.m stack_helm.m)
parameters.min_csgoodratio= 1*ones(length(parameters.periods));%[3 3 3 3 5 10 15 15 15 15 20]; %[3 3 3 3 5 10 15 15]; % minimum radio between good and bad measurements for a good event
parameters.min_phv_tol = 3;
parameters.max_phv_tol = 5;
parameters.is_raydense_weight = 0; %1; % manual says suggested turned off for large azimuthal anisotropy
parameters.min_event_num = 3; %10;
parameters.err_std_tol = 2;
parameters.issmoothmap = 0;%1;
parameters.smooth_wavelength = 0.25;
parameters.event_bias_tol = 3; %2;


%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% parameters for azimuthal anisotropy inversion
% beta version
parameters.smsize = 10; %1;  % averaging nearby grid number
parameters.off_azi_tol = 30; % differ from great circle path in degrees
parameters.is_one_phi = 0;

if length(parameters.periods)~=length(parameters.smweight_array) || length(parameters.periods)~=length(parameters.min_csgoodratio)
    error('Length of periods doesn''t match smweight_array and/or min_csgoodratio');
end

%system(['cp ./setup_parameters.m ',parameters.workingdir]);
