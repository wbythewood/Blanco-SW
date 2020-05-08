% script to take ATaCR code output and put it in the format to be read in
% by ASWMS

% wbh apr 2020

clear;
setup_parameters;

ghDir = '~/Research/github/';
ASWMSDir = [ghDir,'ASWMS-ani/'];
workingDir = parameters.workingdir;
outPath = [workingDir,'eventmat/'];
comp = parameters.component;


ATaCRDir = [ghDir,'ATaCR/'];
evFileName = 'eventtimes_x9testX.txt';
outFileName = [ASWMSDir,'datacache/'];

% get station names
BXHDir = dir(BXHDirName);
xidx = [BXHDir(:).isdir];
xList = {BXHDir(xidx).name};
xList = xList(~ismember(xList,{'.','..'}));

BDHDir = dir(BDHDirName);
didx = [BDHDir(:).isdir];
dList = {BDHDir(didx).name};
dList = dList(~ismember(dList,{'.','..'}));


% get event names
% file name
fid = fopen([ATaCRDir,'config_files/',evFileName]);
% use textscan to read info
evList = textscan(fid,'%s');
%this is for some reason a 1x1 array with 1x[nEvts] cell inside...
evList = evList{1};

% now loop through events
for ie = 1:length(evList)
    % set up matfile name for output
    matFileName = [outPath,char(evList(ie)),'_',comp,'.mat'];
    clear event
    
    % Loop through BXH stations
    % Where are files located
    for ista = 1:length(xList)
        staName = char(xList(ista));
        ifn = [BXHDirName,staName,'/',staName,'_',char(evList(ie)),'_corrseis.mat'];
        disp(ifn)
        staMatFile = dir(fullfile(ifn))
        sta = load(fullfile(staMatFile.folder,staMatFile.name))
        % UGH the ATaCR files don't have event info in them... maybe josh's
        % codes will work better for this after all.
    end
end





