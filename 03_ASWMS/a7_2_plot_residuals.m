%% This script will plot residuals as a function of the azimuth of the incoming ray

% William Hawley
% feb 2021

clear;
setup_parameters
label = 'ALL';
dmin = 60;
minCohere = 0.9;
avgOcePhv = [3.9 3.9 3.9 4.0 4.0 4.0 4.0 4.0 4.1 4.1 4.2];

% paths to various files
workingdir = parameters.workingdir;
matFileDir = parameters.MatFilesDir;
fig_dir_base = parameters.figdir;

% get stations to use
staFileName = parameters.PACStaFile;

% we need dt from the cs files
csFileNames = dir([matFileDir,'CSmeasure/*.mat']);
% and we need isotropic velocities and event info from eikonal step...
%eikonalFileNames = dir([matFileDir,'eikonal_',label,'/*.mat']);

% make a list of the station names we want to use
[nw,staNames,~,~,~] = textread(staFileName,'%s %s %s %s %s');

% the structure we're trying to build should have for each period band, 
% a list of azimuths and delays (% of velocity) so we will need to get 
% azimuth, delay, distance, and average v for each observation. 


dvs = [];
azis = [];

AvgDvs = [];
AvgAzis = [];

% loop over events
for ie = 1:length(csFileNames)
    DvsEvt = [];
    AzisEvt = [];
    temp = load([csFileNames(ie).folder, '/', csFileNames(ie).name]);
    % get station name indices
    staids = [];
    for ista = 1:length(temp.eventcs.stnms)
        if ismember(temp.eventcs.stnms(ista),staNames)
            staids = [staids,ista];
        end
    end
    % get event location
    evLoc = [temp.eventcs.evla,temp.eventcs.evlo];
    % get avg phase vel
    % loop over every observation
    for ics = 1:length(temp.eventcs.CS)
        cs = temp.eventcs.CS(ics);
        % don't look at events if both stations are not in the station list
        if ~ismember(cs.sta1,staids)
            continue
        end
        if ~ismember(cs.sta2,staids)
            continue
        end
        
        % skip if stations too close
        if abs(cs.ddist) < dmin
            continue
        end
        
        % get station locations, use the midpoint as the location to
        % calculate aximuth
        lat1 = temp.eventcs.stlas(cs.sta1);
        lon1 = temp.eventcs.stlos(cs.sta1);
        lat2 = temp.eventcs.stlas(cs.sta2);
        lon2 = temp.eventcs.stlos(cs.sta2);
        [olat,olon] = LatLonMidpoint(lat1,lon1,lat2,lon2);
        [staDist,staAzi] = distance(lat1,lon1,lat2,lon2);
        [dist,azi] = distance(evLoc(1),evLoc(2),olat,olon);
        
        % remove station pairs whose azimuth is too far away from
        % event-array azimuth
        deg_tol = 10; % in degrees
        usedeg = 0;
        if abs(azi - staAzi) < deg_tol
            usedeg = 1;
        elseif abs(abs(azi - staAzi) - 180) < deg_tol
            usedeg = 1;
        elseif abs(abs(azi - staAzi) - 360) < deg_tol
            usedeg = 1;
        end
        if usedeg == 0
            %disp(['bad:'])
            %disp([azi,staAzi,abs(azi-staAzi)])
            continue
            %usedeg = 1;
        end
        %disp(['good'])
        %disp([azi,staAzi,abs(azi-staAzi)])
        
        % now need to loop over periods 
        for ip = 1:length(cs.dtp)
            if ip == 1 && ics == 1
                DvsEvt = zeros(length(temp.eventcs.CS),length(cs.dtp));
                AzisEvt = zeros(length(temp.eventcs.CS),length(cs.dtp));
            end
            % get average phase vel
            avgV = temp.eventcs.avgphv(ip);
            %avgV = avgOcePhv(ip);
            
            % skip if observation is not good
            if cs.isgood(ip) < 1
                continue
            end
            if cs.cohere(ip) < minCohere
                continue
            end
            % need to calculate velocity between the stations from time and
            % distance between them
            vel = cs.ddist / cs.dtp(ip);
            dv = ((vel / avgV) - 1) * 100;
            if dv > 100
                disp(['wow, dv is ',num2str(dv),'%... maybe check this one?'])
                printstr = ['id: ',char(temp.eventcs.id),', Stas: ',char(temp.eventcs.stnms(cs.sta1)),' and ',char(temp.eventcs.stnms(cs.sta2))];
                disp(printstr)
                disp(['ddist = ',num2str(cs.ddist),', dt = ',num2str(cs.dtp(ip)),', ip = ',num2str(ip)])
                continue
            end
            dvs = [dvs,dv];
            azis = [azis,azi];
            Observations(ie).Per(ip).Azi(ics) = azi;
            Observations(ie).Per(ip).Dv(ics) = dv;


        end
    end
    %AvgDvs(ie,:) = DvsEvt;
    %AvgAzis(ie,:) = AzisEvt;
    
end

AvgAzi = [];
AvgDv = [];
for ie = 1:length(Observations)
    for ip = 1:length(Observations(ie).Per)
        AvgAzi(ie,ip) = mean(nonzeros(Observations(ie).Per(ip).Azi));
        AvgDv(ie,ip) = mean(nonzeros(Observations(ie).Per(ip).Dv));
    end
end

figure(71)
clf
hold on
scatter(azis,dvs,'k')

for ip = 1:length(AvgDv(1,:))
    scatter(AvgAzi(:,ip),AvgDv(:,ip),100, ones(1,length(AvgAzi(:,ip)))*ip,'filled')
end
cb = colorbar;
xlabel('Azimuth (\circ)')
ylabel('Phase velocity anomaly (%)')
title(['Anomaly of phase velocity observations, dmin = ',num2str(dmin),' km'])
%ylim([-10 10]);

