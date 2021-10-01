% Apply amplitude correction on the result of eikonal tomography
% also use the stacking result for more accurate correction
% written by Ge Jin, jinwar@gmail.com
% 2013-03

clear;

isfigure = 1;
isoverwrite = 1;

% setup parameters
setup_parameters

% input path and files
% eventcs_path = './CSmeasure/';
% eikonal_data_path = './eikonal/';
% eikonal_stack_file = ['eikonal_stack_',parameters.component];
% helmholtz_path = './helmholtz/';
workingdir = parameters.workingdir;
matFileDir = parameters.MatFilesDir;
eventcs_path = [matFileDir,'CSmeasure/'];
%eikonal_data_path = [matFileDir,'eikonal-flat/'];
eikonal_data_path = [matFileDir,'eikonal/'];
eikonal_stack_file = [matFileDir,'eikonal_stack_',parameters.component];
helmholtz_path = [matFileDir,'helmholtz/'];

figDir = parameters.figdir;
AzFigDir = [figDir,'BackAzHelm/'];

if ~exist(helmholtz_path,'dir')
	mkdir(helmholtz_path);
end

if ~exist(AzFigDir)
    mkdir(AzFigDir);
end

% wbh load plate boundaries
% plate boundaries
mapsDir = [parameters.MapsDir,'PlateBoundaries_NnrMRVL/'];
usgsFN = [parameters.MapsDir,'usgs_plates.txt.gmtdat'];
[pbLat,pbLon] = importPlates(usgsFN);


% load stacked phase velocity map
load(eikonal_stack_file);

% set up useful variables
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode = lalim(1):gridsize:lalim(2);
ynode = lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
amp_var_tol = parameters.amp_var_tol;
alpha_range = parameters.alpha_range;
alpha_search_grid = parameters.alpha_search_grid;
periods = parameters.periods;

eventfiles = dir([eikonal_data_path,'/*_eikonal_',parameters.component,'.mat']);

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
	end
end

for ie = 1:length(eventfiles)
%for ie = 59
	% read in data for this event
	clear eventphv eventcs helmholtz;
	load(fullfile(eikonal_data_path,eventfiles(ie).name));
	eventid = eventphv(1).id;
	matfilename = fullfile(helmholtz_path,[eventphv(1).id,'_helmholtz_',parameters.component,'.mat']);
	if exist(matfilename,'file') && ~isoverwrite
		disp(['exist: ',matfilename,', skip!'])
		continue;
	end
	disp(eventid);
	eventcsfile = [eventcs_path,'/',eventid,'_cs_',parameters.component,'.mat'];
	if exist(eventcsfile,'file')
		load(eventcsfile);
	else
		disp(['Cannot find CS file for ',eventid,', Skipped']);
		continue;
	end
	if length(eventphv) ~= length(eventcs.avgphv)
		disp('Inconsist of period number for CS file and eikonal file');
		continue;
	end

	for ip = 1:length(eventphv)
		%% fit the amplitude surface
		% reset the arrays
		clear stlas stlos amps
		stlas = eventcs.stlas;
		stlos = eventcs.stlos;
		stnms = eventcs.stnms;
		if exist('badstnms','var')
			list_badstaids = find(ismember(eventcs.stnms,badstnms));
		else
			list_badstaids = [];
		end
		amps = zeros(1,length(stlas));
		for ista = 1:length(eventcs.autocor)
			if eventcs.autocor(ista).exitflag(ip)>0
				amps(ista) = eventcs.autocor(ista).amp(ip);
			else
				amps(ista) = NaN;
			end
		end
		% change from power spectrum to amplitude
		amps = amps.^.5;

		% get rid of bad stations
		badstaids = find(isnan(amps));
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		amps(badstaids) = [];
		badstanum = 0; badstaids = [];
		for ista = 1:length(amps)
			if stlas(ista) < lalim(1) || stlas(ista) > lalim(2) || ...
					stlos(ista) < lolim(1) || stlos(ista) > lolim(2) || ismember(ista,list_badstaids);
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
				continue;
			end
			dist = distance(stlas(ista),stlos(ista),stlas,stlos);
			dist = deg2km(dist);
			nearstaids = find(dist > parameters.minstadist & dist < parameters.maxstadist );
			nearstaids(find(ismember(nearstaids,badstaids))) = [];
			if isempty(nearstaids)
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
				continue;
			end
			meanamp = median(amps(nearstaids));
			if amps(ista) < meanamp./amp_var_tol | amps(ista) > meanamp.*amp_var_tol
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
			end
		end
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		amps(badstaids) = [];
		[ampmap,mesh_xi,mesh_yi]=gridfit_jg(stlas,stlos,amps,xnode,ynode,...
							'smooth',2,'regularizer','del4','solver','normal');

		%% Calculate the correction term
		dAmp=del2m(mesh_xi,mesh_yi,ampmap);
		amp_term=-dAmp./ampmap./(2*pi/periods(ip)).^2;
		% smooth the correction term 
		smD=max([300 periods(ip).*parameters.refv]);
		amp_term = gridfit_jg(mesh_xi(:),mesh_yi(:),amp_term(:),xnode,ynode,...
							'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal');
		% prepare the avg phase velocity and event phase velocity
		avgGV = avgphv(ip).GV;
		if sum(size(avgGV)==size(xi)) < 2
			avgGV = interp2(avgphv(ip).xi,avgphv(ip).yi,avgphv(ip).GV,xi,yi,'linear',NaN);
		end
		eventGV = eventphv(ip).GV;
		if sum(size(eventGV)==size(xi)) < 2
			eventxnode = eventphv(ip).lalim(1):eventphv(ip).gridsize:eventphv(ip).lalim(2);
			eventynode = eventphv(ip).lolim(1):eventphv(ip).gridsize:eventphv(ip).lolim(2);
			[eventxi eventyi] = ndgrid(eventxnode,eventynode);
			eventGV = interp2(eventxi,eventyi,eventphv(ip).GV,xi,yi,'linear',NaN);
		end
		% remove the region with small amplitude area
%		meanamp = nanmean(ampmap(:));
%		Tampmap = ampmap';
%		ind = find(Tampmap(:)<meanamp.*parameters.min_amp_tol);
%		eventGV(ind) = NaN;
		% apply correction
		[GV_cor alpha_errs alphas] = amp_correct(avgGV, eventGV, amp_term, alpha_range, alpha_search_grid);
        
        % JBR - line to replace negative phase velocities (imaginary) with event velocity
        GV_cor(imag(GV_cor)~=0) = eventGV(imag(GV_cor)~=0);
        
		[temp bestalphai] = min(alpha_errs);
		bestalpha = alphas(bestalphai);
		fprintf('%f ',bestalpha);

		% fill in informations
		helmholtz(ip).evla = eventphv(ip).evla;
		helmholtz(ip).evlo = eventphv(ip).evlo;
		helmholtz(ip).raydense = eventphv(ip).raydense;
		helmholtz(ip).goodnum = eventphv(ip).goodnum;
		helmholtz(ip).badnum = eventphv(ip).badnum;
		helmholtz(ip).id = eventphv(ip).id;
		helmholtz(ip).xi = xi;
		helmholtz(ip).yi = yi;
		helmholtz(ip).GV_cor = GV_cor;
		helmholtz(ip).GV = eventGV;
		helmholtz(ip).alpha_errs = alpha_errs;
		helmholtz(ip).alphas = alphas;
		helmholtz(ip).bestalpha = bestalpha;
		helmholtz(ip).amp_term = amp_term;
		helmholtz(ip).ampmap = ampmap;
		helmholtz(ip).period = periods(ip);
		bestalphas(ip,ie) = bestalpha;

		% plot to check
		if isfigure
			figure(37)
			clf
                        set(gcf,'renderer','zbuffer');
			subplot(2,2,1)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,eventGV);
            %load seiscmap %wbh
			%colormap(gca,seiscmap) %wbh
			colormap(gca,roma) %wbh
            cb1 = colorbar; %wbh
            ylabel(cb1,'Vs (km/s)') %wbh
			%colorbar
			title('before cor');
			subplot(2,2,2)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,GV_cor);
			cb2 = colorbar; %wbh
			title('after cor');
			nanind = find(isnan(eventGV(:)));
			ampmap = ampmap';
			ampmap(nanind) = NaN;
			amp_term = amp_term';
			amp_term(nanind) = NaN;
			subplot(2,2,3)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,ampmap);
			title('amplitude map')
			plotm(stlas,stlos,'v')
            colormap default %wbh
			cb3 = colorbar; %wbh
			subplot(2,2,4)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,amp_term);
			cb4 = colorbar; %wbh
			[temp bestalphai] = min(alpha_errs);
			title('correction term')
            drawnow;
            
			figure(38)
			clf
                        set(gcf,'renderer','zbuffer');
			plot(alphas,alpha_errs,'x');
            drawnow;
%             pause;
%             % wbh figure for all periods like a6
%             N=3; 
%             M = floor(length(periods)/N) +1;
%             Mphvel = floor((1+length(periods))/N) +1;   % wbh
%             figure(39)
%             clf
%             for ip = 1:length(periods)
%                 subplot(M,N,ip)
%                 ax = worldmap(lalim,lolim);
%                 h1 = surfm(xi,yi,GV_cor);
%                 title(["Period: ",num2str(periods(ip))],'fontsize',15)
%                 avgv = nanmean(eventphv(ip).GV(:));
%                 if isnan(avgv)
%                     continue;
%                 end
%                 r = 0.1;
%                 caxis([avgv*(1-r) avgv*(1+r)])
%                 colorbar
%                 load seiscmap
%                 colormap(seiscmap)
%                 h = colorbar; %wbh
%                 ylabel(h,'Vs (km/s)') %wbh
%             end
		end % end of isfigure
	end  % loop of period
    % wbh figure for all periods like a6
    N=3; 
    M = floor(length(periods)/N) +1;
    Mphvel = floor((1+length(periods))/N) +1;   % wbh
    figure(39)
    clf
    for ip = 1:length(periods)
        %subplot(M,N,ip)
        subplot(Mphvel,N,ip)
        ax = worldmap(lalim,lolim);
        h1 = surfm(xi,yi,GV_cor);
        title(["Period: ",num2str(periods(ip))],'fontsize',15)
        %avgv = nanmean(eventphv(ip).GV(:));
        avgv = nanmean(helmholtz(ip).GV_cor(:));
        if isnan(avgv)
            continue;
        end
        r = 0.2;
        caxis([avgv*(1-r) avgv*(1+r)])
        colorbar
        %load seiscmap
        %colormap(seiscmap)
        colormap(roma)
        h = colorbar; %wbh
        ylabel(h,'Vs (km/s)') %wbh

    end
    % wbh add station map
    subplot(Mphvel,N,length(periods)+1)
        
        ax = worldmap(lalim,lolim);
        set(ax, 'Visible', 'off')
        plotm(pbLat,pbLon,'LineWidth',2,'Color','k') % plate boundaries
        %geoshow(eventcs.stlas,eventcs.stlos,'DisplayType','point','Marker','.','MarkerEdgeColor','b','MarkerSize',12)
        amps = [];
        for ista = 1:length(eventcs.stnms)
            amps(ista) = mean(eventcs.autocor(ista).amp);
        end
        scatterm(eventcs.stlas,eventcs.stlos,30,amps,'o','filled')
        caxis([min(amps) max(amps)]);
        h = colorbar;
        ylabel(h,'Amplitude');
        % get event azimuth
        az = azimuth(eventcs.stlas(1),eventcs.stlos(1),eventcs.evla,eventcs.evlo);
        [distdeg,az] = distance(eventcs.stlas(1),eventcs.stlos(1),eventcs.evla,eventcs.evlo);
        arrlen = 1.4;
        arru = arrlen.*cosd(az);
        arrv = arrlen.*sind(az);
        arrLat = mean(lalim);
        arrLon = mean(lolim);
        quiverm(arrLat,arrLon,arru,arrv,'r')
        azstr = ["az: "+num2str(round(az))+'\circ'];
        textm(lalim(1)+0.5,lolim(1)+0.5,azstr,'FontSize',12) 
    %sgtitle("Corrected Phase Velocities for "+eventcs.id);
    
    MwStr = sprintf('%.2f',eventphv(ip).Mw);
    DistStr = sprintf('%.0f',distdeg);
    sgtitle("Corrected Phase Velocities for "+eventcs.id+' M'+MwStr+' Dist: '+DistStr+'\circ'); %wbh

    %ofn = [figDir,eventid,'/PhaseVels_0.25_HCorr.png'];   % wbh save in event dir
    ofn = [AzFigDir,num2str(round(az)),'_',DistStr,'_M',MwStr,'_HelmPhaseVels_0.25.png'];   % wbh save in BackAz dir
    %set(gcf,'PaperUnits','centimeters');
    %set(gcf,'PaperPosition', [0,0,10,15]);
    saveas(gcf,ofn) 
    drawnow;
	matfilename = fullfile(helmholtz_path,[eventphv(1).id,'_helmholtz_',parameters.component,'.mat']);
	save(matfilename,'helmholtz');
	fprintf('\n');
	disp(['Saved to ',matfilename]);
end  % loop of events
