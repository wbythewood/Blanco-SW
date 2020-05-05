% make_starttimes
%
% Make start time files for daily data downloads from list of event times in
% format YYYYmmddhhMM
%
% J. Russell & H. Janiszewski
% hjaniszewski@carnegiescience.edu
% updated 11/19

% W.B. Hawley update to use parameter file
% updated 02/20

setup_parameter;

% Load event list
evlist = textread(evFile,'%s');

fid = fopen(dayFile,'w');
for iev = 1:length(evlist)
    evnum = datenum(evlist(iev),'yyyymmddHHMM');
    daynums = flip(evnum - [1:Ndays]);
    for iday = 1:length(daynums)
        day = datestr(daynums(iday),'yyyymmddHHMM');
        fprintf(fid,'%s\n',day);
    end
end
fclose(fid);
