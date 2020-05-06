% correctevent.m
% script for correcting event files for tilt and compliance noise

clear all; close all
setup_parameter;

isfigure_sta = 1;
isfigure_spectra = 1;
issavefigure = 1;
isoverwrite =1;

T1 = 10; T2= 150; %filter period range for plotting seismic data

inpath_event = EventPreproDir; %'../../data/datacache_prepro/'; % path to earthquake data (data to be corrected)

% input lists of stations with bad components, file doesn't have to exists
% leave as default if no bad stations
badstalist = '../../Bad_Z.txt';
badhorzlist = '../../Bad_H.txt';
badpreslist = '../../Bad_P.txt';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
setup_parameter;

if exist(badstalist)
	badstas = textread(badstalist,'%s');
else
    badstas = [];
end
if exist(badhorzlist)
	badhorz = textread(badhorzlist,'%s');
else
    badhorz = [];
end
if exist(badpreslist)
	badpres = textread(badpreslist,'%s');
else
    badpres = [];
end

for itf = 1:length(TF_list)
    disp([itf])
    TF_size(itf) = length(TF_list{itf});
end
[vec,TF_list_sort] = sort(TF_size,'descend');

event_filename = dir(fullfile(inpath_event));

inpath_dir = sprintf('%s/TRANSFUN/',OUTdir);
TF_stadir = dir(fullfile(inpath_dir));
ii=1;
for itf = 1:length(TF_stadir)
    station = TF_stadir(itf).name;
    inpath = sprintf('%s/TRANSFUN/%s/',OUTdir,station);
    if ~isdir(inpath)
        continue
    end
    filenames = dir(fullfile(inpath,['*.mat']));
    if isempty(filenames)
        continue
    end
    TF_stalist{ii} = station;
    ii=ii+1;
end

for ie = 1:length(event_filename)
    close all
    clear corrseis_matrix
    inpath = sprintf('%s%s/',inpath_event,event_filename(ie).name);
    if ~isdir(inpath)
        continue
    end
    if length(event_filename(ie).name)~=12;
        continue
    end
    eventid = event_filename(ie).name;
    disp(['Event: ',eventid]);
    event_info_dir = dir(fullfile(inpath_event,event_filename(ie).name,['*.mat']));
    stas = length(event_info_dir);

    for ista = 1:stas %loop through stations for given event
        clear spec_mat
        TF_check = zeros(size(TF_size));
        idx1 = findstr(event_info_dir(ista).name,'_');
        idx2 = findstr(event_info_dir(ista).name,'.mat');

        network = event_info_dir(ista).name(idx1(1)+1:idx1(2)-1);
        station = event_info_dir(ista).name(idx1(2)+1:idx2-1)

        netsta = [network,station];

        TF_staidx = find(strcmp(netsta,TF_stalist)==1);
        if isempty(TF_staidx)
            continue
        end

        disp(netsta);

        sta = load(fullfile(inpath_event,event_filename(ie).name,'/',event_info_dir(ista).name)); %loading the correct file
        disp([inpath_event,event_filename(ie).name,'/',event_info_dir(ista).name]); %display file name for check, can turn off

        % Specify transfer function file path and output paths
        inpath_trans = sprintf('%s/TRANSFUN/%s/',OUTdir,netsta);
        if tf_op ==1
            outpath = sprintf('%s/CORRSEIS/%s/',OUTdir,netsta);
            figoutpath=sprintf('%s/CORREVENTS',FIGdir);
        elseif tf_op ==2
            outpath = sprintf('%s/CORRSEISAVTF/%s/',OUTdir,netsta);
            figoutpath=sprintf('%s/CORREVENTS',FIGdir);
        end

		if ~isoverwrite && exist(sprintf('%s/%s_%s_corrseis.mat',outpath,netsta,eventid))==2
            disp([eventid,' exists. Skipping...'])
            continue
        end

        [Zraw,H1raw,H2raw,Praw,taxisZ,taxis1,taxis2,taxisP,dt] = varsetup_correctevent(sta,chz_vec,ch1_vec,ch2_vec,chp_vec,T);
        if isnan(Zraw);
            continue
        end


        % FINDING AND LOADING THE TF FILES
        trans_filename=findTFs_correctevent(inpath_trans,tf_op,netsta,eventid);
        if isnan(trans_filename)
            continue
        end

        if exist(trans_filename,'file')
            goodP = 1;
            goodH=1;
            goodZ=1;
            if ~isempty(badstas)
            if ~isempty(find(~cellfun('isempty', strfind(badstas,netsta)))) == 1
                disp('Vertical Bad.')
                goodZ=0; %if pressure is marked as bad, turn good flag off
            end
            end
            if ~isempty(badpres)
            if ~isempty(find(~cellfun('isempty', strfind(badpres,netsta)))) == 1
                disp('Pressure Bad.')
                goodP=0; %if pressure is marked as bad, turn good flag off
            end
            end
            if ~isempty(badhorz)
            if ~isempty(find(~cellfun('isempty', strfind(badhorz,netsta)))) == 1
                disp('Horizontal Bad.')
                goodH=0; %if horizontal is marked as bad, turn good flag off
            end
            end
        end

        % Calculate event spectra properties
        [spec_mat,npad0,npts,f,NFFT]=spectra_correctevent(Zraw,H1raw,H2raw,Praw,dt,taxisZ,taxis1,taxis2,taxisP,taptime);

        disp('Calculating Corrected Seismograms...')
        disp(sprintf('Using TF: %s', trans_filename));

        if ~exist(figoutpath)
            mkdir(figoutpath);
        end
        if ~exist(outpath)
            mkdir(outpath);
        end

        %%%%%%%%%%%% LOAD TRANSFER FUNCTIONS
        load(trans_filename);

        freqcomp=transprop.params.freqcomp;
        if filop==1
            lp = taper_lim(2);
            hp = taper_lim(1);
        elseif filop ==2
            lp = freqcomp+0.005;
            hp=0;
        end

        NFFT = transprop.params.NFFT;
        f=transprop.params.f;
        for itf = 1:length(TFs); %taper each of the individual transfer functions
            TFs(itf).transfunc_tap = tf_taper((TFs(itf).transfunc)',f',hp,lp,tap_width);
        end



        ii = 1;
        for itf = 1:length(TF_list_sort)
            TF_cal = TF_list(TF_list_sort(itf)); % name of the TF we're applying now
            if TF_check(TF_list_sort(itf)) == 0 % check if TF calculated already
                % check if rotational or not, i.e., H max
                if ~isempty(strfind(cell2mat(TF_cal),'H')) % if rotational
                    if  tf_op == 2 % average doesn't make sense for this
                        continue
                    end
                    if length(TF_cal{1}) == 2 % 1 component rotational TF, i.e. ZH
                        [corrspec,corrtime,corrgood,label_list] = tfcomp1rotate_correctevent(TF_cal,TFs,transprop,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                        TF_check(TF_list_sort(itf)) = 1;
                        corrseis(ii).label = label_list{1};
                        corrseis(ii).spectrum = corrspec;
                        corrseis(ii).timeseries = corrtime;
                        corrseis(ii).isgood = corrgood;
                        ii = ii+1;
                    elseif length(TF_cal{1}) == 4 % ZP-H
                        [corrspec1_2,corrspec1_32,corrtime1_2,corrtime1_32,corrgood1_2,corrgood1_32,...
                            label_list] = tfcomp2rotate_correctevent(TF_cal,TFs,transprop,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                        TF_check(TF_list_sort(itf)) = 1;
                        corrseis(ii).label = label_list{1};
                        corrseis(ii).spectrum = corrspec1_2;
                        corrseis(ii).timeseries = corrtime1_2;
                        corrseis(ii).isgood = corrgood1_2;
                        if ~isempty(find(strcmp(label_list{1},TF_list)==1))
                            TFidx = find(strcmp(label_list{1},TF_list)==1);
                            TF_check(TFidx)=1;
                        end
                        corrseis(ii+1).label = label_list{2};
                        corrseis(ii+1).spectrum = corrspec1_32;
                        corrseis(ii+1).timeseries = corrtime1_32;
                        corrseis(ii+1).isgood = corrgood1_32;
                        if ~isempty(find(strcmp(label_list{2},TF_list)==1))
                            TFidx = find(strcmp(label_list{2},TF_list)==1);
                            TF_check(TFidx)=1;
                        end
                        ii = ii+2;
                    end

                 else % if component wise, i.e., H1 and H2
                    if length(TF_cal{1}) == 2 % 1 component correction
                        [corrspec,corrtime,corrgood,label_list] = tfcomp1_correctevent(TF_cal,TFs,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                        TF_check(TF_list_sort(itf)) = 1;
                        corrseis(ii).label = label_list{1};
                        corrseis(ii).spectrum = corrspec;
                        corrseis(ii).timeseries = corrtime;
                        corrseis(ii).isgood = corrgood;
                        ii = ii+1;
                    elseif length(TF_cal{1}) == 5 % 3 component correction, i.e., ZP-21
                        [corrspec1_2,corrspec1_32,corrspec1_432,corrtime1_2,...
                            corrtime1_32,corrtime1_432,corrgood1_2,corrgood1_32,corrgood1_432,label_list] =...
                            tfcomp3_correctevent(TF_cal,TFs,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                        TF_check(TF_list_sort(itf)) = 1;
                        corrseis(ii).label = label_list{1}; % this is Z1 correction...
                        corrseis(ii).spectrum = corrspec1_2;
                        corrseis(ii).timeseries = corrtime1_2;
                        corrseis(ii).isgood = corrgood1_2;
                        if ~isempty(find(strcmp(label_list{1},TF_list)==1))
                            TFidx = find(strcmp(label_list{1},TF_list)==1);
                            TF_check(TFidx)=1;
                        end
                        corrseis(ii+1).label = label_list{2}; % Z2-1 correction
                        corrseis(ii+1).spectrum = corrspec1_32;
                        corrseis(ii+1).timeseries = corrtime1_32;
                        corrseis(ii+1).isgood = corrgood1_32;
                        if ~isempty(find(strcmp(label_list{2},TF_list)==1))
                            TFidx = find(strcmp(label_list{2},TF_list)==1);
                            TF_check(TFidx)=1;
                        end
                        corrseis(ii+2).label = label_list{3}; % and ZP-21 correction
                        corrseis(ii+2).spectrum = corrspec1_432;
                        corrseis(ii+2).timeseries = corrtime1_432;
                        corrseis(ii+2).isgood = corrgood1_432;
                        if ~isempty(find(strcmp(label_list{3},TF_list)==1))
                            TFidx = find(strcmp(label_list{3},TF_list)==1);
                            TF_check(TFidx)=1;
                        end
                        ii=ii+3;
                    end
                end
            end
        end

        % wbh edits?try applying instrument response here
        % first need pz from the z component raw data
        % '/Users/whawley/Research/github/ATaCR/data/datacache/201201231604'
        staFileName = sprintf('%s_%s_%s.mat',eventid,network,station);
        %load(fullfile(EventDataDir,eventid,staFileName))
        load(fullfile(EventPreproDir,eventid,staFileName))
%         dataDirName = dir(fullfile(EventDataDir));
%         %dateFileName = sprintf('%s%s/',EventDataDir,event_filename(fileLoop).name);
%         dateFileName = sprintf('%s%s/',EventPreprocDir,event_filename(fileLoop).name);
%         stationFileNames = dir(fullfile(inpath,['*.mat']));
%         if isempty(stationFileNames)
%             continue
%         end
%         for ista = 1:length(stationFileNames)
%             if isempty(findstr([NetworkName,'_',station],stationFileNames(ista).name))
%                 continue
%             end
%             load(fullfile(EventDataDir,'/',dataDirName(fileLoop).name,'/',stationFileNames(ista).name));
%
%         end
        zIndex = find(contains({traces.channel},'Z'));
        zTrace = traces(zIndex);
        for compi = 1:length(corrseis)
            trace.data = corrseis(compi).timeseries;
            trace.network = network;
            trace.station = station;
            trace.channel = zTrace.channel;
            trace.sacpz = zTrace.sacpz;
            trace.sampleCount = length(corrseis(compi).timeseries);
            trace.sampleRate = trace.sampleCount/EventDataLength;
            trace.sensitivity = zTrace.sensitivity;
            trace.sensitivityFrequency = zTrace.sensitivityFrequency;
            trace.sensitivityUnits = zTrace.sensitivityUnits;
            if RespAfterFlag == 1
                NewTrace = rm_resp(trace,eventid,LoPassCorner,nPoles,PZDir,1,length(corrseis));
                corrseis(compi).timeseries = NewTrace.data_cor;
            end
        end

        corrected.params.f = f;
        corrected.params.tf_op = tf_op;
        corrected.params.eventid = eventid;
        corrected.params.filop = filop;
        corrected.params.taxis = taxisZ;
        corrected.params.station = station;
        corrected.params.network = network;
        corrected.params.NFFT = NFFT;
        corrected.params.dt = dt;
        corrected.params.TFfilename = trans_filename;
%         corrected.params.stalat =
%         corrected.params.stalon =
%         corrected.params.freqcomp =

        filename = sprintf('%s/%s_%s_corrseis.mat',outpath,netsta,eventid);
        save(filename,'corrseis','corrected');

        %station plots
        if isfigure_sta == 1;
            spitplots_correctevent(dt,T1,T2,Zraw,H1raw,H2raw,Praw,taxisZ,eventid,netsta,corrseis,f,NFFT)
            if issavefigure==1
                figure(101)
                filename=sprintf('%s/%s_%s_originalseis',figoutpath,eventid,netsta);
                print(gcf,'-dpng',filename)
                if tf_op ==1
                    figure(102)
                    filename=sprintf('%s/%s_%s_corrseis',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                elseif tf_op==2
                    figure(102)
                    filename=sprintf('%s/%s_%s_corrseis_av',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                end
                figure(103)
                filename=sprintf('%s/%s_%s_corrspectra',figoutpath,eventid,netsta);
                print(gcf,'-dpng',filename)

                figure(104)
                filename=sprintf('%s/%s_%s_seis+spectra',figoutpath,eventid,netsta);
                print(gcf,'-dpng',filename)
            end
        end
		%plot spectra against Peterson's low and high noise models
		if isfigure_spectra
            plot_corrspectra(T1,T2,Zraw,corrseis,f,NFFT,dt,taxisZ,eventid,netsta,lp,hp)
            if issavefigure==1
                if tf_op ==1
                    figure(1)
                    filename=sprintf('%s/%s_%s_spectra',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                elseif tf_op==2
                    figure(1)
                    filename=sprintf('%s/%s_%s_spectra_av',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                end
            end
        end
    end

end
