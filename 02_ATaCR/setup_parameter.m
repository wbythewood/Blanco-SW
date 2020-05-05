% parameters for the NoiseTiltComp_pck

% path for matlab codes and functions
addpath ('function');
javaaddpath('IRIS-WS-2.0.18.jar');

%--- Directory Structure ---%

%Base Directory for Output
BaseDir = 'X9_test2D/';%''X9_M6.5';
% location of unpreprocessed matlab files
NoiseDataDir = strcat(BaseDir,'data/noise_day/'); % output folder for data
% location of the continuous matlab files for spectral properties
NoisePreproDir = strcat(BaseDir,'data/noise_day_prepro/');
% location of unpreprocessed event data
EventDataDir = strcat(BaseDir,'data/event_data/');
% location of processed event data
EventPreproDir = strcat(BaseDir,'data/event_data_prepro/');
% location of pole zero directory
PZDir = ''; % leave blank if no PZs
%PZDir = strcat(BaseDir,'data/PZDir');
% output directory for spectra
OUTdir = strcat(BaseDir,'data/NOISETC/');
% directory for figure output
FIGdir = strcat(BaseDir,'figures/');

% paths for the event and noise time lists
evFile = 'config_files/eventtimes_X9test2.txt';
dayFile = 'config_files/starttimes_X9test2.txt';


%evFile = 'config_files/Events_X9_M6.5.txt';
%dayFile = 'config_files/starttimes_X9_M6.5.txt.txt';

%--- Data to download ---%

% networks -- still only can do one
%NetworkNames = '7D';
NetworkName = 'X9';
% stations
%StationNames = {'M07A'};
%StationNames = {'*'};
StationNames = textread('config_files/X9_D_stations.txt','%s');

% Response Removal
% option of removing response from Z component only after corrections have 
% been applied -- 0 is do not remove response after, 1 is remove response 
% after correcting

RespAfterFlag = 1;

% Pre-processing high-pass filter
% this will apply a high pass filter to the raw data
% useful to remove very long period signals from the data
FilterBeforeFlag = 1;
LoPassCorner = 0.005; %Hz

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
EventDataLength = 7200;

%--- Preprocessing ---%

% Sample Rate of seismic data you want -- will downsample to this
%SampleRate = 5;
SampleRate = 1; %LHZ is 1 sps
% preprocessing high-pass filter (if hpFilt above is 1)
nPoles = 5;

%--- Downloading Specrtral Properties ---%

% Spectral Properties Windowing
% the legnth of each time window, in sec, should match the event data 
% length for corrections
T    = 7200; 

% fraction of window overlap for spectra calculation
overlap = 0.3; 

% Quality Control Parameters for Daily Windows
pb = [0.004 .2]; % pass-band, in Hz
tolerance = 1.5; % tolerance factor for QC
a_val = 0.05;    % f-test for QC (1 - a_val = confidence)
minwin = 10;     % minimum numbers of time window for segment to be accepted

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
