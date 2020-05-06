% download_event

% downloads the event data files that will be corrected for tilt and
% compliance noise

% H. Janiszewski
% hjaniszewski@carnegiescience.edu
% updated 11/17
% W.B. Hawley updated to include paramter file
% updated 02/20

clear; close all
setup_parameter;

datacache = EventDataDir; % output folder for data

%%%%% end user input parameters %%%%%

if ~exist(datacache,'dir')
    mkdir(datacache)
end

startlist = textread(evFile,'%s');
chanlist = [chz_vec, ch1_vec, ch2_vec, chp_vec];
chanlist = strjoin(chanlist,',');

for id = 1:length(startlist)
   eventid = cell2mat(startlist(id));
   disp(sprintf('Start Time: %s',eventid));
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
       error = 0;
       stnm = stations_info(ista).StationCode;
       network = stations_info(ista).NetworkCode;
       if ~exist(fullfile(datacache,eventid),'dir')
           mkdir(fullfile(datacache,eventid));
       end
       sta_filename = fullfile(datacache,eventid,[eventid,'_',network,'_',stnm,'.mat']);
       if exist(sta_filename,'file')
           disp(['Exist: ',sta_filename,', Skip!']);
           continue;
       end
       disp(['Downloading station: ',stnm,' From:',starttime,' To:',endtime]);
		try
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
