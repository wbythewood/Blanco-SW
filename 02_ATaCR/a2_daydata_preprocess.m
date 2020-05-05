% data_preprocess
% Collection of sample data pre-processing codes that may be applied to
% data prior to noise analysis and tilt/compliance correction. This script
% will apply the pre-processing steps and save the data in an appropriate
% file structure.
% Pre-procssing steps included here are:
%  - Downsampling
%  - Response Removal
%  - Gain Adjustment
% Each may be included or excluded depending on the needs of the user.
% HAJ July 2016

% W.B. Hawley updated to include paramter file
% also updated to do event and noise calculations in the same code
% updated 02/20

clear; close all
setup_parameter;

% check if directories exist
if ~exist(NoisePreproDir)
    mkdir(NoisePreproDir);
end
if ~exist(EventPreproDir)
    mkdir(EventPreproDir);
end

% common parameters
network = NetworkName;
stations = StationNames;
channels = [chz_vec, ch1_vec, ch2_vec, chp_vec];

% flags
RemoveResp = [chz_resp, ch1_resp, ch2_resp, chp_resp];
GainCorr = [chz_gain, ch1_gain, ch2_gain, chp_gain];
samprate = SampleRate;
HiPassFilt = [chz_hpFilt, ch1_hpFilt, ch2_hpFilt, chp_hpFilt];

% now that have network, make dir inside prepro dir
if ~exist(fullfile(NoisePreproDir,network),'dir')
    mkdir(fullfile(NoisePreproDir,network));
end
if ~exist(fullfile(EventPreproDir,network),'dir')
    mkdir(fullfile(EventPreproDir,network));
end

%% Preprocess Noise Data
for ista = 1:length(stations)
    station = stations{ista};
    if ~exist(fullfile(NoisePreproDir,network,station),'dir')
        mkdir(fullfile(NoisePreproDir,network,station));
    end
    % data_filenames = dir(fullfile(EventDataDir,network,station,'/',['*.mat']));
    filesuff = sprintf('*_%s_%s.mat',network,station);
    data_filenames = dir(fullfile(NoiseDataDir,network,station,'/',filesuff));
    for ie =1 : length(data_filenames) % begin file loop for each station
        load(fullfile(NoiseDataDir,network,station,'/',data_filenames(ie).name));
        traces_day_new = traces_day;
        dayid = data_filenames(ie).name(1:12);
        prob=0;

        for ic = 1:length(channels) % begin channel loop
            chan = channels(ic);
            idxch = find(ismember({traces_day.channel},chan));
            if length(idxch)>1
                disp('Skipping. Too many records for single channel.')
                prob = 1;
                continue
            end
            if isempty(idxch)
                continue
            end
            chan_data = traces_day(idxch);
            rate = chan_data.sampleRate;
            dt = 1/rate;
            nsp = 2;
            %wbh add a high-pass filter
            if FilterBeforeFlag == 1
                fNyquist = 0.5 * chan_data.sampleRate;
                FreqHiPass = 3600; %(s)
                [B,a] = butter(2,1/FreqHiPass/fNyquist,'high');
                chan_data.data = filtfilt(B,a,chan_data.data);
            end
            %%%%%%%%%%%%%%%%%
            % REMOVE RESPONSE
            %%%%%%%%%%%%%%%%%
            if RemoveResp(ic) ==1
                if isempty(traces_day(1).sacpz.poles) & isempty(PZDir)
                    prob = 1;
                    continue
                end
                chan_data = rm_resp(chan_data,dayid,LoPassCorner,nPoles,PZDir,idxch,nsp);
                data_raw = chan_data.data_cor;
            else
                data_raw = chan_data.data;
            end
            %%%%%%%%%%%%%%%%%
            % GAIN CORRECTION
            %%%%%%%%%%%%%%%%%
            data_raw = GainCorr(ic).*data_raw;
            %%%%%%%%%%%%%%%%%
            % HIGH PASS FILTER - should link this and the bit in rm_resp
            % function to a different function, so they are the same, but the
            % user can specify parameters easily
            %%%%%%%%%%%%%%%%%
            if HiPassFilt(ic) == 1

                lo_w=2*pi*LoPassCorner;

                N = length(data_raw);
                delta = dt;
                Tr = N*delta;

                if mod(N,2)
                    faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/Tr);
                else
                    faxis = [0:N/2,-N/2+1:-1]*(1/Tr);
                end
                w = faxis.*2*pi;

                hpfiltfrq=( ((w./lo_w).^(2*nPoles))./(1+(w./lo_w).^(2*nPoles)) );
                norm_trans=hpfiltfrq;    % this is normalization transfer function
                norm_trans(find(isnan(norm_trans))) = 0;

                fftdata = fft(data_raw);
                fftdata = fftdata(:).*norm_trans(:);
                data_raw = real(ifft(fftdata));
            end
            %%%%%%%%%%%%%%%%%
            % DOWNSAMPLING
            %%%%%%%%%%%%%%%%%
            if rate == samprate
                taxis = [0:dt:(length(data_raw)-1)*dt]';
            else
                dt_new = 1/samprate;
                [data_raw,taxis] = resample(data_raw,dt_new,dt);
            end

            traces_day_new(idxch).data = data_raw;
            traces_day_new(idxch).sampleRate = samprate;
            traces_day_new(idxch).sampleCount = length(data_raw);
        end % end channel loop
        %save good files
        if prob ==0
            traces_day = traces_day_new;
            filename = fullfile(NoisePreproDir,network,station,'/',data_filenames(ie).name);
            save(filename,'traces_day');
        end
    end % end file loop
end % end station loop

%% Preprocess Event Data

data_filenames = dir(fullfile(EventDataDir));
for ista = 1:length(stations)
    station = stations{ista};
    for ie = 1 : length(data_filenames) % begin file loop for each station
        inpath = sprintf('%s%s/',EventDataDir,data_filenames(ie).name);
        if ~isdir(inpath)
            continue
        end
        if length(data_filenames(ie).name)~=12;
            continue
        end
        station_filenames = dir(fullfile(inpath,['*.mat']));
        if isempty(station_filenames)
            continue
        end
        for is = 1:length(station_filenames)
            if isempty(findstr([network,'_',station],station_filenames(is).name))
                continue
            end
            load(fullfile(EventDataDir,'/',data_filenames(ie).name,'/',station_filenames(is).name));
        traces_new = traces;
        eventid = data_filenames(ie).name(1:12);
        if ~exist([EventPreproDir,eventid])
        mkdir([EventPreproDir,eventid]);
    end
        prob=0;
        for ic = 1:length(channels) % begin channel loop
            chan = channels(ic);
            idxch = find(ismember({traces.channel},chan));
            if length(idxch)>1
                disp('Skipping. Too many records for single channel.')
                prob = 1;
                continue
            end
            if isempty(idxch)
                continue
            end
            chan_data = traces(idxch);
            rate = chan_data.sampleRate;
            dt = 1/rate;
            % debug strange taper behavior
            nsp = 2;
             if idxch == 4
                 subplot(nsp,1,1)
                 t = (chan_data.sampleCount/chan_data.sampleRate);
                 time = linspace(0,t,chan_data.sampleCount);
                 data = chan_data.data;
                 plot(time,data,'k');
                 xlabel('Time (s)');
                 title('Z raw');
                 xlim([min(time),max(time)]);
             end

            % try to do a really long-period high-pass to get rid of super
            % long-period noise?
            if FilterBeforeFlag == 1
                fNyquist = 0.5 * chan_data.sampleRate;
                FreqHiPass = 3600; %(s)
                [B,a] = butter(2,1/FreqHiPass/fNyquist,'high');
                chan_data.data = filtfilt(B,a,chan_data.data);
            end
%             if idxch == 4
%                  subplot(nsp,1,2)
%                  data = chan_data.data;
%                  plot(time,data,'k');
%                  xlabel('Time (s)');
%                  title('Z hi-pass 3600s');
%                  %title('same as above');
%                  xlim([min(time),max(time)]);
%              end
            %%%%%%%%%%%%%%%%%
            % REMOVE RESPONSE
            %%%%%%%%%%%%%%%%%
            if RemoveResp(ic) ==1
                %chan_data = rm_resp(chan_data,eventid,LoPassCorner,nPoles,PZDir);
                chan_data = rm_resp(chan_data,eventid,LoPassCorner,nPoles,PZDir,idxch,nsp);
                data_raw = chan_data.data_cor;
            else
                data_raw = chan_data.data;
            end
%             if idxch == 4
%                  subplot(4,1,2)
%                  plot(time,data_raw,'k');
%                  xlabel('Time (s)');
%                  title('Z Resp Removed');
%                  xlim([min(time),max(time)]);
%              end
            %%%%%%%%%%%%%%%%%
            % GAIN CORRECTION
            %%%%%%%%%%%%%%%%%
            data_raw = GainCorr(ic).*data_raw;
%             if idxch == 4
%                  subplot(4,1,3)
%                  plot(time,data_raw,'k');
%                  xlabel('Time (s)');
%                  title('Z Gain Corr');
%                  xlim([min(time),max(time)]);
%              end
             %%%%%%%%%%%%%%%%%
            % HIGH PASS FILTER - should link this and the bit in rm_resp
            % function to a different function, so they are the same, but the
            % user can specify parameters easily
            %%%%%%%%%%%%%%%%%
            if HiPassFilt(ic) ==1

            lo_w=2*pi*LoPassCorner;

            N = length(data_raw);
            delta = dt;
            Tr = N*delta;

            if mod(N,2)
                faxis = [0:(N-1)/2,-(N-1)/2:-1]*(1/Tr);
            else
                faxis = [0:N/2,-N/2+1:-1]*(1/Tr);
            end
            w = faxis.*2*pi;

            hpfiltfrq=( ((w./lo_w).^(2*nPoles))./(1+(w./lo_w).^(2*nPoles)) );
            norm_trans=hpfiltfrq;    % this is normalization transfer function
            norm_trans(find(isnan(norm_trans))) = 0;

            fftdata = fft(data_raw);
            fftdata = fftdata(:).*norm_trans(:);
            data_raw = real(ifft(fftdata));
            end
            if idxch == 4
                 subplot(nsp,1,nsp)
                 %t = (chan_data.sampleCount/chan_data.sampleRate);
                 %time = linspace(0,t,chan_data.sampleCount);
                 %data_prepro = chan_data.data;
                 plot(time,data_raw,'k');
                 xlabel('Time (s)');
                 title('Z Preprocessed');
                 xlim([min(time),max(time)]);
             end
            %%%%%%%%%%%%%%%%%
            % DOWNSAMPLING
            %%%%%%%%%%%%%%%%%
            if rate == samprate
                taxis = [0:dt:(length(data_raw)-1)*dt]';
            else
                dt_new = 1/samprate;
                [data_raw,taxis] = resample(data_raw,dt_new,dt);
            end

            traces_new(idxch).data = data_raw;
            traces_new(idxch).sampleRate = samprate;
            traces_new(idxch).sampleCount = length(data_raw);
            % debug strange taper behavior

        end % end channel loop
        %save good files
        if prob ==0
            traces = traces_new;
            filename = fullfile(EventPreproDir,eventid,'/',station_filenames(is).name);
            save(filename,'traces');
        end
        end
    end % end file loop
end % end station loop
