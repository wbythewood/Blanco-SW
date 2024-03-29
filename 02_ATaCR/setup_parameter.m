% parameters for the Blanco-SW/02_ATaCR

% path for matlab codes and functions
addpath ('function');
javaaddpath('./IRIS-WS-2.0.18.jar');

%--- Directory Structure ---%

%Base Directory for Output
% old mac
BaseDir = '/Users/wbhawley/Research/Seismology/Blanco-SW/';
% gaherty mac
BaseDir = '/Users/whawley/Research/github/Blanco-SW/';
% new mac
BaseDir = '/Users/whawley/Research/Blanco-SW/';

DataDir = strcat(BaseDir,'data/');
AtacrProcDir = strcat(BaseDir,'02_ATaCR/dataProc/');
%Label = 'GordaAddition_test';
Label = 'CoverAzimuth';

% If Downloading using MATLAB code...
% location of unpreprocessed matlab files
%NoiseDataDir = strcat(BaseDir,'data/noise_day/'); % output folder for data
% location of the continuous matlab files for spectral properties
%NoisePreproDir = strcat(BaseDir,'data/noise_day_prepro/');
%NoisePreproDir = strcat(BaseDir,'data/noise_day_prepro_',Label,'/');
% location of unpreprocessed event data
%EventDataDir = strcat(BaseDir,'data/event_data/');
% location of processed event data
%EventPreproDir = strcat(BaseDir,'data/event_data_prepro/');
%EventPreproDir = strcat(BaseDir,'data/noise_day_prepro_',Label,'/');

% wbh editing to keep intermediate processing step files separate
NoiseDataDir = strcat(AtacrProcDir,'data/',Label,'/noise_day/');
NoisePreproDir = strcat(AtacrProcDir,'data/',Label,'/noise_day_prepro/');
EventDataDir = strcat(AtacrProcDir,'data/',Label,'/event_data/');
EventPreproDir = strcat(AtacrProcDir,'data/',Label,'/event_data_prepro/');
% location of pole zero directory
PZDir = ''; % leave blank if no PZs
%PZDir = strcat(BaseDir,'data/PZDir');

% If already downloaded using SAC...
% location of SAC noise files
sacDayData = strcat(DataDir,'SAC_Noise_',Label,'/');
% location of SAC event files
sacEventData = strcat(DataDir,'SAC_Events_',Label,'/');

% output directory for spectra
OUTdir = strcat(BaseDir,'data/NOISETC/',Label,'/');
% directory for figure output
FIGdir = strcat(BaseDir,'figures/ATaCR/',Label,'/');

% paths for the event and noise time lists
%evFile = strcat(BaseDir,'config/BlancoProblematicEvt.txt');
%dayFile = strcat(BaseDir,'config/BlancoProblematicDay.txt');
% for removing from ambient noise... test 05 jan 
% same file for "events" and noise... since no events
evFile = strcat(BaseDir,'config/BlancoEvents_',Label,'_M6.5.txt');
dayFile = strcat(BaseDir,'config/BlancoNoise_',Label,'_M6.5.txt');

evFile = strcat(BaseDir,'config/BlancoEvents_',Label,'_M6.3-5.0.txt');
dayFile = strcat(BaseDir,'config/BlancoNoise_',Label,'_M6.3-5.0.txt');

%--- Data to download ---%

% networks -- still only can do one
% stations
%StationNames = {'*'};

%NetworkName = 'X9';
%StationNames = textread(strcat(BaseDir,'config/X9_stations.txt'),'%s');
%StationNames = textread(strcat(BaseDir,'config/X9_stations_ANTtest.txt'),'%s');

NetworkName = '7D';
%StationNames = textread(strcat(BaseDir,'config/7D_stations.txt'),'%s');
StationNames = textread(strcat(BaseDir,'config/7D_stations_withGorda.txt'),'%s');
%StationNames = textread(strcat(BaseDir,'config/7D_stations_ANTtest.txt'),'%s');

% Response Removal
% option of removing response from Z component only after corrections have
% been applied -- 0 is do not remove response after, 1 is remove response
% after correcting

RespAfterFlag = 0;

% Pre-processing high-pass filter
% this will apply a high pass filter to the raw data
% useful to remove very long period signals from the data
FilterBeforeFlag = 1;
LoPassCorner = 0.005; %Hz
LoPassCorner = 0.001; %Hz
LoPassCorner = 1/30;

% Channel Names and Corrections
% in the a_ sections there can only be one vertical channel per station
% downloaded... maybe I will go back to fix this later. When do we want
% more than one downloaded... especially if in the b_ sections, it skips if
% there's more than one component? Just need the flexibility.

% *_vec is the channel names
% *_resp is whether to remove response for above channels - 0 = no; 1 = yes
% *_gain is gain correction - multiply data by this number
% *_hpFilt will apply a high-pass filter -- should only be 1 for _resp=0

% Channel Info
chz_vec = {'LHZ'};
chz_resp = [1];
chz_gain = [1];
chz_hpFilt = [0];

ch1_vec = {'LH1'};
ch1_resp = [1];
ch1_gain = [1];
ch1_hpFilt = [0];

ch2_vec = {'LH2'};
ch2_resp = [1];
ch2_gain = [1];
ch2_hpFilt = [0];

%chp_vec = {'BXH'};
chp_vec = {'BDH'};
chp_resp = [1];
chp_gain = [1];
chp_hpFilt = [0];

% if removing response after, don't need to remove resp before...
if RespAfterFlag == 0
    chz_resp = zeros(size(chz_resp));
    ch1_resp = zeros(size(ch1_resp));
    ch2_resp = zeros(size(ch2_resp));
    chp_resp = zeros(size(chp_resp));
end

% Timing Info
% number of noise days to use for calculating spectra
Ndays = 4;
% length of record (seconds) for data download
NoiseDataLength = 86400;
EventDataLength = 6000;
% ANT no events... 
%EventDataLength = 86400;
%EventDataLength = 43200;
%NoiseDataLength = EventDataLength;

%--- SAC I/O ---%
% which corrected seismogram to use:
% 1=Z1; 2=Z2-1; 3=ZP-21; 4=ZH; 5=ZP-H; 6=ZP
corr_idx = 3;

%--- Preprocessing ---%

% Sample Rate of seismic data you want -- will downsample to this
SampleRate = 1; %LHZ is 1 sps
% preprocessing high-pass filter (if hpFilt above is 1)
nPoles = 5;

%--- Downloading Specrtral Properties ---%

% Spectral Properties Windowing
% the legnth of each time window, in sec, should match the event data
% length for corrections

T = 'Use EventDataLength'; %unclear why these weren't the same before...

% fraction of window overlap for spectra calculation
overlap = 0.3;
%ant
%overlap = 0.01;

% Quality Control Parameters for Daily Windows
pb = [0.004 .2]; % pass-band, in Hz
tolerance = 1.5; % tolerance factor for QC
a_val = 0.05;    % f-test for QC (1 - a_val = confidence)
minwin = 10;     % minimum numbers of time window for segment to be accepted
%minwin = 1;

% Tilt orientation - only matters if using transfer functions with the 'H'
% option, but package needs variables specified to run; leave as default if
% not using
tiltfreq = [.005, .035]; % specifying frequency ranges for maximum coherence search
% tiltfreq = [1/20, 1/5]; % specifying frequency ranges for maximum coherence search

% Quality Control Parameters for Deployment Days
pb_dep = [0.004 .2]; %pass-band, in Hz
tolerance_dep = 1.5;
a_val_dep = 0.05;

% Transfer Function Options
TF_list = {'ZP','Z1','Z2-1','ZP-21','ZH','ZP-H'};
% convention: Z1 = transfer function between Z and H1
%             Z2-1 = transfer function between Z with H1 removed, and H2
%		      ZP-21 = transfer function between Z with H1 and H2 removed, and P
%             ZH = transfer function between Z and rotated max horizontal noise direction (H)

% Correction Options
taptime = 0.075; % taper for seismogram correction step

tf_op = 1; %option for using either daily (1) or station average (2) TF in correction

filop = 2; %how to filter TF
% 1 - user specified constant high pass and low pass
% 2 - %lowpass - 0.005+freqcomp, adaptive to the infragravity wave cutoff, no high pass;
tap_width = 0.01; %width in percent of frequency vector of cosine taper is this actually used?????
taper_lim = [0 1]; % upper and lower freuqncy in Hz limit for taper if user specified (option 1); zero means not applied


% %% Octave Things
% pkg load financial
% pkg load signal