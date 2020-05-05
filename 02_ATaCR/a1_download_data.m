% download_data

% downloads the data files used in calculating the noise spectra and the
% transfer functions for tilt and compliance corrections and saves them as
% matfiles (default is 24 hours of data in each file).

% H. Janiszewski
% hjaniszewski@carnegiescience.edu
% updated 2/18

% W.B. Hawley updated to include paramter file
% also download event data in the same script
% updated 02/20

clear; close all
setup_parameter;

% check if directories exist
if ~exist(NoiseDataDir,'dir')
    mkdir(NoiseDataDir)
end
if ~exist(EventDataDir,'dir')
    mkdir(EventDataDir)
end

% a few parameters that are similar across noise+event data
chanlist = [chz_vec, ch1_vec, ch2_vec, chp_vec];
chanlist = strjoin(chanlist,',');

%% first download the noise data

startlist = textread(dayFile,'%s');


for id = 1:length(startlist)
   eventid = cell2mat(startlist(id));
   fprintf('Start Time: %s\n',eventid);
   otime = datenum(eventid,'yyyymmddHHMM');
   starttime = datestr(otime,'yyyy-mm-dd HH:MM:SS');
   endtime = datestr(otime+NoiseDataLength/3600/24,'yyyy-mm-dd HH:MM:SS');

   stations_info = irisFetch.Stations('channel',NetworkName,StationNames,'*',chz_vec,'startTime',starttime,'endTime',endtime);
   disp([stations_info])

   for ista =1:length(stations_info)
       error = 0;
       disp([stations_info(ista)])
       stnm = stations_info(ista).StationCode;
       network = stations_info(ista).NetworkCode;
       if ~exist(fullfile(NoiseDataDir,network),'dir')
           mkdir(fullfile(NoiseDataDir,network));
       end
       if ~exist(fullfile(NoiseDataDir,network,stnm),'dir')
           mkdir(fullfile(NoiseDataDir,network,stnm));
       end
       sta_filename = fullfile(NoiseDataDir,network,stnm,[eventid,'_',network,'_',stnm,'.mat']);
%        if exist(sta_filename,'file')
%            disp(['Exist: ',sta_filename,', Skip!']);
%            continue;
%        end
       disp(['Downloading station: ',stnm,' From:',starttime,' To:',endtime]);
		   try
            traces_day = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime,'includePZ');
			 save(sta_filename,'traces_day');
		   catch e
            e.message;
            error = 1;
       end
        if error ==1
            try % to try and get around the missing zeros for some pressure components
                %traces_day = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime);
                %save(sta_filename,'traces_day');
            catch e
                e.message;
                continue;
            end
        end
    end

end
%% Now the event data

startlist = textread(evFile,'%s');

for id = 1:length(startlist)
   eventid = cell2mat(startlist(id));
   fprintf('Start Time: %s\n',eventid);
   if length(eventid) == 12
        otime = datenum(eventid,'yyyymmddHHMM'); % look at this when I get back
   elseif length(eventid) == 14
       otime = datenum(eventid,'yyyymmddHHMMSS');
       eventid = eventid(1:12); % for naming purposes only, start time will still be saved to the second in traces file
   end
   starttime = datestr(otime,'yyyy-mm-dd HH:MM:SS');
   endtime = datestr(otime+EventDataLength/3600/24,'yyyy-mm-dd HH:MM:SS');

   stations_info = irisFetch.Stations('channel',NetworkName,StationNames,'*',chz_vec,'startTime',starttime,'endTime',endtime);


   for ista =1:length(stations_info)
       %disp([stations_info(ista)])
       error = 0;
       stnm = stations_info(ista).StationCode;
       network = stations_info(ista).NetworkCode;
       if ~exist(fullfile(EventDataDir,eventid),'dir')
           mkdir(fullfile(EventDataDir,eventid));
       end
       sta_filename = fullfile(EventDataDir,eventid,[eventid,'_',network,'_',stnm,'.mat']);
%        if exist(sta_filename,'file')
%            disp(['Exist: ',sta_filename,', Skip!']);
%            continue;
%        end
       disp(['Downloading station: ',stnm,' From:',starttime,' To:',endtime]);
		try
            % disp([network,stnm,chanlist,starttime,endtime])
            traces = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime,'includePZ');
            save(sta_filename,'traces');
		catch e
            e.message;
            error = 1;
        end
        if error ==1
            try % to try and get around the missing zeros for some pressure components
                traces = irisFetch.Traces(network,stnm,'*',chanlist,starttime,endtime);
                save(sta_filename,'traces');
            catch e
                e.message;
                continue;
            end
        end
    end

end
