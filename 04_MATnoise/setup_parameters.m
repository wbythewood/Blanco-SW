% Setup_parameters for MATnoise
%
% NJA, 4/2/2016
% JBR, 6/16/2016
% WBH, 11/15/2020

addpath('./functions/');
addpath('./functions/calc_Rayleigh_disp/');

IDString = '10-30s_mingrv2';

%%% --- Paths to important files --- %%%
parameters.workingdir = '/Users/whawley/Research/github/Blanco-SW/';
parameters.NoiseDir = [parameters.workingdir,'/04_MATnoise/'];
parameters.dataDir = [parameters.workingdir,'data/'];
parameters.configDir = [parameters.workingdir,'config/'];
% path to data that has NOT yet been corrected for t/c noise
parameters.datapath = [parameters.dataDir,'Sac_Noise_Test/'];
% and path to data that HAS been corrected for t/c noise
%parameters.datapath = [parameters.dataDir,'CORRSEIS_SAC/'];
parameters.PZpath = [parameters.dataDir,'PZ/'];
parameters.StaListFile = [parameters.configDir,'stalist_nw.txt'];
%parameters.ccfpath = [parameters.NoiseDir,'CCF/',IDString,'/'];
parameters.ccfpath = [parameters.NoiseDir,'CCF/Full_raw/'];
parameters.figpath = [parameters.workingdir,'figures/ANT/',IDString,'/'];  % wbh addition
parameters.seis_path = [parameters.NoiseDir,'seismograms/'];
parameters.orientation_path = [parameters.dataDir,'orientations.txt'];
%[stalist, stalat, stalon, staz] = textread(parameters.StaListFile,'%s %f %f %f\n');
[nwlist, stalist, stalat, stalon, staz] = textread(parameters.StaListFile,'%s %s %f %f %f\n'); % wbh add nw parameter
parameters.nwlist = nwlist;
parameters.stalist = stalist;
parameters.stalat = stalat;
parameters.stalon = stalon;
parameters.staz = staz;
parameters.nsta = length(parameters.stalist);

%%% --- Parameters to build up gaussian filters --- %%%
parameters.min_width = 0.18;
parameters.max_width = 0.30;

%%% --- Parameters for initial processing --- %%%
parameters.dt = 1; % sample rate
parameters.comp = 'LH'; % component
parameters.mindist = 20; % min. distance in kilometers

%%% --- Parameters for a1_ccf_ambnoise --- %%%
parameters.strSACcomp = 'Z';
parameters.strNAMEcomp = 'ZZ';
parameters.winlength = 3; %hours
parameters.winDirName = ['window',num2str(parameters.winlength),'hr'];
parameters.Nstart_sec = 50; % number of sections to offset start of seismogram
parameters.IsRemoveIR = 0; % remove instrument response
parameters.units_RemoveIR = 'M'; % 'M' displacement | 'M/S' velocity
parameters.IsDetrend = 1; % detrend the data
parameters.IsTaper = 1; % Apply cosine taper to data chunks
% (1) ONE-BIT NORMALIZATION & SPECTRAL WHITENING? (Bensen et al. 2007)
parameters.IsSpecWhiten = 0; % Whiten spectrum
parameters.IsOBN = 1; % One-bit normalization
% (2) TIME-FREQUENCY NORMALIZATION (Ekstrom et al. 2009; Shen et al. 2011)
parameters.IsFTN = 0; % Frequency-time normalization? (If 1, applied instead of whitening and one-bit normalization)
parameters.frange_FTN = [1/60 1/5]; % frequency range over which to construct FTN seismograms
% (3) BASIC PREFILTER (Ekstrom 2011)
parameters.IsPrefilter = 0; % apply butterworth bandpass filter before cross-correlation?
parameters.frange_prefilt = [1/50 1/3]; % note in FREQ (1/T)

%%% --- Parameters for a2_plot_ccf_record --- %%%
parameters.PeriodRange = [10 30]; % note in PERIOD
parameters.trace_space = 0; % km
parameters.snr_thresh = 2.5;
parameters.dep_tol = [0 0]; % [sta1, sta2] OBS Depth tolerance;
parameters.max_grv = inf; %5.5;
parameters.min_grv = 2.0; %1.4; %1.6
parameters.xlims = [-250 250];
parameters.ylims = [0 450];
parameters.h20_grv = 1.5;
parameters.costap_wid = 0.2; % 0 => box filter; 1 => Hann window
parameters.IsButterworth = 0; % 1=butterworth; 0=tukey
if parameters.IsButterworth == 1
    parameters.FiltStr = 'butterworth';
elseif parameters.IsButterworth == 0
    parameters.FiltStr = 'tukey';
else
    disp('Error with parameters.IsButterworth');
end

%%% --- Parameters for a6_fitbessel --- %%%
parameters.npts = parameters.winlength*3600;
parameters.Wavelengths = 0; % number of wavelengths to fit wbh 12/2020 appears to break if nonzero...
parameters.npers = 8;    % number of periods
parameters.damp = [1; 1; 1]; % [fit, smoothness, slope]
parameters.is_normbessel = 0; % normalize bessel function by analytic envelope?
parameters.iswin = 1; % Use the time-domain windowed ccfs?

%%% --- Parameters for using Radon Transform picks --- %%%
parameters.path_LRT_picks = './mat-LRTdisp/LRT_picks/';


%%% --- Parameters for b_ Tomography --- %%%
addpath('./tomo_functions');
parameters.lalim = [41.5 45.5];%[42.5 45] ;
parameters.lolim = [-131.5 -125.25];%[-131.5 -125.5];
parameters.gridsize = 0.5; %0.25;   % in degrees
parameters.gridsize_azi = 0.5; %3; %1.5; % gridsize for 2D azimuthal anisotropy (degrees)

% Smoothing parameters
parameters.smweight0 = 100; % isotropic second derivative smoothing
parameters.smweight0_azi = 1e3; %1000; % anisotropic second derivative smoothing
parameters.flweight0_azi = 1000; %1000; % anisotropic first derivative flatness

parameters.raydensetol=deg2km(parameters.gridsize)*0.25; %deg2km(parameters.gridsize); %deg2km(parameters.gridsize)*2;
parameters.raydensetol_azi=deg2km(parameters.gridsize_azi)*0.25; %deg2km(parameters.gridsize)*2;
parameters.fiterrtol = 2;   % error allowed in the wavelet fitting
parameters.dterrtol = 4;    % largest variance of the inversion error allowed
parameters.maxerrweight = 5; % Maximum error weight
parameters.polyfit_dt_err = 2; % (s) dt error greater than this, weighted 0
parameters.tomo_snr_tol = 1.5; %2.5; % minimum signal-to-noise
parameters.r_tol_min = 20; % [km] minimum station separation
parameters.r_tol_max = 600; % [km] maximum station separation
parameters.err_tol = 1.5; %0.5; % maximum misfit of bessel fit between observed and synthetic

parameters.is_rtolmin_wavelength = 0; % determine distance tolerance by wavelength?
parameters.wl_fac = 1.0; % wavelength of above

parameters.r = 0.03; %0.01; % controls color bar [avgv(1-r) avgv(1+r)]
