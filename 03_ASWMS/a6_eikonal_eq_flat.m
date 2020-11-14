% Read in the eventcs structures and apply eikonal tomography on each event.
% Written by Ge Jin, jinwar@gmail.com
% 2013.1.16
%
% JBR 3/19/19 : add first derivative smoothing 
%
clear
%plot native

% setup parameters
setup_parameters
setup_ErrorCode

% JBR
cohere_tol = parameters.cohere_tol;

fiterr_tol = 1e-2; % wavelet fitting error, throw out measurements greater than this
maxstadist = 600;
minstadist = 200;
maxstadist = parameters.maxstadist;
minstadist = parameters.minstadist;
%cohere_tol = 0.80;
%min_stadist_wavelength = 0.33; %0.5; % minimum station separation in wavelengths
min_stadist_wavelength = 0; % minimum station separation in wavelengths
max_stadist_wavelength = 999;
ref_phv = [3.9936 4.0041 4.0005 3.9999 3.9929 3.9832 3.9813 3.9841 3.9874 3.9996 4.0138 4.0519 4.0930 4.1677 4.2520]; % for calculating wavelength

% debug setting
isfigure = 1;
isdisp = 0;
is_overwrite = 1;

% % input path
% eventcs_path = './CSmeasure/';
% % output path
% eikonl_output_path = './eikonal/';

workingdir = parameters.ASWMSDir;
matFileDir = parameters.MatFilesDir;
% input path
eventcs_path = [matFileDir,'CSmeasure/'];
% output path
eikonl_output_path = [matFileDir,'eikonal-flat/'];
% figures directory
fig_dir_base = parameters.figdir;
figDirPhv = [parameters.figdir,'BackAz-flat/'];
badStaList = [parameters.configDir,'badsta.lst'];


if ~exist(figDirPhv)
    mkdir(figDirPhv);
end

if ~exist(eikonl_output_path)
	mkdir(eikonl_output_path);
end

% plate boundaries
mapsDir = [parameters.MapsDir,'PlateBoundaries_NnrMRVL/'];
usgsFN = [parameters.MapsDir,'usgs_plates.txt.gmtdat'];
[pbLat,pbLon] = importPlates(usgsFN);

comp = parameters.component;
lalim=parameters.lalim;
lolim=parameters.lolim;
gridsize=parameters.gridsize;
periods = parameters.periods;
raydensetol=parameters.raydensetol;
smweight_array = parameters.smweight_array;
flweight_array = parameters.flweight_array; % JBR
Tdumpweight0 = parameters.Tdumpweight;
Rdumpweight0 = parameters.Rdumpweight;
fiterrtol = parameters.fiterrtol;
dterrtol = parameters.dterrtol;
isRsmooth = parameters.isRsmooth;
inverse_err_tol = parameters.inverse_err_tol;
min_amp_tol  = parameters.min_amp_tol;

% setup useful variables
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

% Setup universal smoothing kernel
disp('initial the smoothing kernel')
tic
	% longtitude smoothing
    [i,j] = ndgrid(1:Nx,2:(Ny-1));
    ind = j(:) + Ny*(i(:)-1);
    dy = diff(ynode)*cosd(mean(xnode));  % correct smoothing for latitude
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
                    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny,Nx*Ny);

	% latitude smoothing
    [i,j] = ndgrid(2:(Nx-1),1:Ny);
    ind = j(:) + Ny*(i(:)-1);
    dx = diff(xnode);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
            [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny,Nx*Ny)];

    F=sparse(Nx*Ny*2*2,Nx*Ny*2);
    for n=1:size(Areg,1)
        ind=find(Areg(n,:)~=0);
        F(2*n-1,2*ind-1)=Areg(n,ind);
        F(2*n,2*ind)=Areg(n,ind);
    end
toc

% JBR - define first derivative "flatness" kernel
F2 = flat_kernel_build(xnode, ynode, Nx*Ny);

% read in bad station list, if existed
% wbh use config file, dont hardwire
if exist(badStaList)
	badstnms = textread(badStaList,'%s');
	disp('Found Bad stations:')
	disp(badstnms)
end
% if exist('badsta.lst')
% 	badstnms = textread('badsta.lst','%s');
% 	disp('Found Bad stations:')
% 	disp(badstnms)
% end

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
for ie = 1:length(csmatfiles)
%for ie = 30
	clear eventphv 
	% read in data and set up useful variables
	temp = load([eventcs_path,csmatfiles(ie).name]);
	eventcs =  temp.eventcs;
	disp(eventcs.id)
	evla = eventcs.evla;
	evlo = eventcs.evlo;
    

    % set up fig dir name
    figDir = [fig_dir_base,eventcs.id,'/'];
    if ~exist(figDir)
        mkdir(figDir);
    end

	matfilename = [eikonl_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat'];
	if exist(matfilename,'file') && ~is_overwrite
		disp(['Exist ',matfilename,', skip']);
		continue;
	end


	if exist('badstnms','var')
		badstaids = find(ismember(eventcs.stnms,badstnms));
	else
		badstaids = [];
	end

	% Build the rotation matrix
	razi = azimuth(xi+gridsize/2,yi+gridsize/2,evla,evlo)+180;
	R = sparse(2*Nx*Ny,2*Nx*Ny);
	for i=1:Nx
		for j=1:Ny
			n=Ny*(i-1)+j;
			theta = razi(i,j);
			R(2*n-1,2*n-1) = cosd(theta);
			R(2*n-1,2*n) = sind(theta);
			R(2*n,2*n-1) = -sind(theta);
			R(2*n,2*n) = cosd(theta);
		end
	end

	% Calculate the relative travel time compare to one reference station
	travel_time = Cal_Relative_dtp_wbh(eventcs,badstaids);

	% Build the ray locations
	clear rays 
	for ics = 1:length(eventcs.CS)
		rays(ics,1) = eventcs.stlas(eventcs.CS(ics).sta1);
		rays(ics,2) = eventcs.stlos(eventcs.CS(ics).sta1);
		rays(ics,3) = eventcs.stlas(eventcs.CS(ics).sta2);
		rays(ics,4) = eventcs.stlos(eventcs.CS(ics).sta2);
	end

	% Build the kernel
	disp('Buildling up ray path kernel')
	tic
		mat=kernel_build(rays,xnode,ynode);
	toc

	% build dumping matrix for ST
	dumpmatT = R(2:2:2*Nx*Ny,:);
	
	% build dumping matrix for SR
	dumpmatR = R(1:2:2*Nx*Ny-1,:);
	
	% Loop through the periods
	for ip = 1:length(periods)
		smweight0 = smweight_array(ip);
        flweight0 = flweight_array(ip); % JBR
		dt = zeros(length(eventcs.CS),1);
		w = zeros(length(eventcs.CS),1);
        ddist = zeros(length(eventcs.CS),1);
        N_coh_tol = 0;
        N_minWavelength = 0;
        N_maxWavelength = 0;
        N_maxDist = 0;
        N_minDist = 0;
        N_fiterr = 0;
        N_usable = 0;
        N_isbadBefore = 0;
        N_badSta = 0;
		for ics = 1:length(eventcs.CS)
            if eventcs.CS(ics).cohere(ip)<cohere_tol && eventcs.CS(ics).isgood(ip)>0
                eventcs.CS(ics).isgood(ip) = ErrorCode.low_cohere;
                N_coh_tol = N_coh_tol + 1;
            end
			%if (eventcs.CS(ics).ddist < ref_phv(ip)*periods(ip)*min_stadist_wavelength || ...
			if (abs(eventcs.CS(ics).ddist) < ref_phv(ip)*periods(ip)*min_stadist_wavelength || ...
                    eventcs.CS(ics).ddist > ref_phv(ip)*periods(ip)*max_stadist_wavelength) && eventcs.CS(ics).isgood(ip)>0
				eventcs.CS(ics).isgood(ip) = ErrorCode.min_stadist_wavelength;
                %N_wavelength = N_wavelength + 1;
                if eventcs.CS(ics).ddist > ref_phv(ip)*periods(ip)*max_stadist_wavelength
                    N_maxWavelength = N_maxWavelength + 1;
                end
                if abs(eventcs.CS(ics).ddist) < ref_phv(ip)*periods(ip)*min_stadist_wavelength
                    N_minWavelength = N_minWavelength + 1;
                    disp([num2str(eventcs.CS(ics).ddist),' < ',num2str(ref_phv(ip)*periods(ip)*min_stadist_wavelength),' (wavelength)'])
                    disp((['Stations ',eventcs.stnms(eventcs.CS(ics).sta1),' ',eventcs.stnms(eventcs.CS(ics).sta2)]))
                    disp(' ')
                end
            end
            if (eventcs.CS(ics).ddist > maxstadist || ...
                    eventcs.CS(ics).ddist < minstadist) && eventcs.CS(ics).isgood(ip)>0
                eventcs.CS(ics).isgood(ip) = -13;
                %N_maxDist = N_maxDist + 1;
                if eventcs.CS(ics).ddist > maxstadist
                    N_maxDist = N_maxDist + 1;
                end
                %if eventcs.CS(ics).ddist < minstadist
                if abs(eventcs.CS(ics).ddist) < minstadist
                    N_minDist = N_minDist + 1;
                    disp([num2str(eventcs.CS(ics).ddist),' < ',num2str(minstadist)])
                    disp(['Stations ',eventcs.stnms(eventcs.CS(ics).sta1),' ',eventcs.stnms(eventcs.CS(ics).sta2)])
                    disp(' ')
                end
            end
            if (eventcs.CS(ics).fiterr(ip) > fiterr_tol)
                eventcs.CS(ics).isgood(ip) = -14;
                N_fiterr = N_fiterr + 1;
            end
            if eventcs.CS(ics).cohere(ip)>=cohere_tol && eventcs.CS(ics).isgood(ip)==ErrorCode.low_cohere
                eventcs.CS(ics).isgood(ip) = 1;
            end
			if eventcs.CS(ics).isgood(ip) > 0 
				dt(ics,:) = eventcs.CS(ics).dtp(ip);
				w(ics,:) = 1;
                N_usable = N_usable + 1;
			else
				dt(ics,:) = eventcs.CS(ics).dtp(ip);
				w(ics,:) = 0;
                N_isbadBefore = N_isbadBefore + 1;
			end
			if sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids)) > 0
				w(ics,:) = 0;
                N_badSta = N_badSta + 1;
            end
            ddist(ics,:) = eventcs.CS(ics).ddist;
        end
        disp(['Too low coherence: ', num2str(N_coh_tol)])
        disp(['Stations too many wavelengths apart: ',num2str(N_maxWavelength)])
        disp(['Stations too few wavelengths apart: ',num2str(N_minWavelength)])
        disp(['Stations too far apart: ',num2str(N_maxDist)])
        disp(['Stations too close: ',num2str(N_minDist)])
        disp(['Fit error too high: ',num2str(N_fiterr)])
        disp(['Number of good observations: ',num2str(N_usable)])
        disp(['Number of prior bad observations: ',num2str(N_isbadBefore)])
        disp(['Number of observations with bad stations: ',num2str(N_badSta)])
		W = sparse(diag(w));

		% Normalize smoothing kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;
        
        % JBR - Normalize flatness kernel
        NR=norm(F2,1);
        NA=norm(W*mat,1);
        flweight = flweight0*NA/NR;

		% Normalize dumping matrix for ST
		NR=norm(dumpmatT,1);
		NA=norm(W*mat,1);
		dumpweightT = Tdumpweight0*NA/NR;
		
		% Normalize dumping matrix for SR
		NR=norm(dumpmatR,1);
		NA=norm(W*mat,1);
		dumpweightR = Rdumpweight0*NA/NR;

		% Set up matrix on both side
		if isRsmooth
            A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
        else
            A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
        end

		avgv = eventcs.avgphv(ip);
        rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
        
		% Least square inversion
        phaseg=(A'*A)\(A'*rhs);
	        
        % Iteratively down weight the measurement with high error
		niter=0;
		ind = find(diag(W)==0);
		if isdisp
			disp(['Before iteration'])
			disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
			disp(['Bad Measurement Number: ', num2str(length(ind))]);
		end
        niter=1;
        while niter < 2
            niter=niter+1;
            err = mat*phaseg - dt;
			err = W*err;
            %            err = W*err;
            stderr=std(err);
            if stderr > dterrtol
                stderr = dterrtol;
            end
            for i=1:length(err)
                if abs(err(i)) > inverse_err_tol*stderr
                    W(i,i)=0;
                end
            end
            ind = find(diag(W)==0);
			if isdisp
				disp('After:')
				disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
				disp(['Bad Measurement Number: ', num2str(length(ind))]);
			end
            
            % Rescale the smooth kernel
            NR=norm(F,1);
            NA=norm(W*mat,1);
            smweight = smweight0*NA/NR;
            
            % JBR - Normalize flatness kernel
            NR=norm(F2,1);
            NA=norm(W*mat,1);
            flweight = flweight0*NA/NR;
            
            % rescale dumping matrix for St
            NR=norm(dumpmatT,1);
            NA=norm(W*mat,1);
            dumpweightT = Tdumpweight0*NA/NR;
            
            % rescale dumping matrix for SR
            NR=norm(dumpmatR,1);
            NA=norm(W*mat,1);
            dumpweightR = Rdumpweight0*NA/NR;
            
            if isRsmooth
                A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
            else
                A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
            end
            rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
            phaseg=(A'*A)\(A'*rhs);
        end	

        % Calculate the kernel density
        %sumG=sum(abs(mat),1);
        ind=1:Nx*Ny;
        rayW = W;
        rayW(find(rayW>1))=1;
        raymat = rayW*mat;
        sumG(ind)=sum((raymat(:,2*ind).^2+raymat(:,2*ind-1).^2).^.5,1);
        clear raydense
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                raydense(i,j)=sumG(n);
            end
        end
        
        %        disp(' Get rid of uncertainty area');
        fullphaseg = phaseg;
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                if raydense(i,j) < raydensetol %&& ~issyntest
                    phaseg(2*n-1)=NaN;
                    phaseg(2*n)=NaN;
                end
            end
        end

		% Change phaseg into phase velocity
		for i=1:Nx
			for j=1:Ny
				n=Ny*(i-1)+j;
				GVx(i,j)= phaseg(2*n-1);
				GVy(i,j)= phaseg(2*n);
			end
		end
		GV=(GVx.^2+GVy.^2).^-.5;
		% Get rid of uncertain area
        % Forward calculate phase velocity
        phv_fwd = ddist./(mat*phaseg(1:Nx*Ny*2));

		% save the result in the structure
		eventphv(ip).rays = rays;
		eventphv(ip).w = diag(W);
		eventphv(ip).goodnum = length(find(eventphv(ip).w>0));
		eventphv(ip).badnum = length(find(eventphv(ip).w==0));
		eventphv(ip).dt = dt;
		eventphv(ip).GV = GV;
		eventphv(ip).GVx = GVx;
		eventphv(ip).GVy = GVy;
        eventphv(ip).phv_fwd = phv_fwd;
		eventphv(ip).raydense = raydense;
		eventphv(ip).lalim = lalim;
		eventphv(ip).lolim = lolim;
		eventphv(ip).gridsize = gridsize;
		eventphv(ip).id = eventcs.id;
		eventphv(ip).evla = eventcs.evla;
		eventphv(ip).evlo = eventcs.evlo;
		eventphv(ip).evdp = eventcs.evdp;
		eventphv(ip).period = periods(ip);
		eventphv(ip).traveltime = travel_time(ip).tp;
		eventphv(ip).stlas = eventcs.stlas;
		eventphv(ip).stlos = eventcs.stlos;
		eventphv(ip).stnms = eventcs.stnms;
        eventphv(ip).isgood = eventphv(ip).w>0;
		eventphv(ip).Mw = eventcs.Mw;
		disp(['Period:',num2str(periods(ip)),', Goodnum:',num2str(eventphv(ip).goodnum),...
				'Badnum:',num2str(eventphv(ip).badnum)]);
	end % end of periods loop
	if isfigure
		N=3; M = floor(length(periods)/N) +1;
        Mphvel = floor((1+length(periods))/N) +1;   % wbh

		figure(88)
		clf
		for ip = 1:length(periods)
			%subplot(M,N,ip)
            subplot(Mphvel,N,ip)   % wbh
			ax = worldmap(lalim, lolim);
			set(ax, 'Visible', 'off')
			h1=surfacem(xi,yi,eventphv(ip).GV);
			% set(h1,'facecolor','interp');
%			load pngcoastline
%			geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)
			title(['Periods: ',num2str(periods(ip))],'fontsize',15)
			avgv = nanmean(eventphv(ip).GV(:));
			if isnan(avgv)
				continue;
			end
			r = 0.1;
			caxis([avgv*(1-r) avgv*(1+r)])
			colorbar
			load seiscmap
			colormap(seiscmap)
            h = colorbar; %wbh
            ylabel(h,'Vs (km/s)') %wbh
            
%             % Plot ray paths
%             jj=0;
%             for ii = 1:length(eventphv(ip).w)
%                 if ~eventphv(ip).isgood(ii) %|| raytomo(ip).sta_dep(ixsp) > -4500
%                     continue
%                 end
%                 jj = jj + 1;
%                 lat1(jj,1) = eventphv(ip).rays(ii,1);
%                 lon1(jj,1) = eventphv(ip).rays(ii,2);
%                 lat2(jj,1) = eventphv(ip).rays(ii,3);
%                 lon2(jj,1) = eventphv(ip).rays(ii,4);
%             end
%             hold on;
%             h = plotm([lat1 lat2]',[lon1 lon2]','-k','linewidth',0.5); hold on;
		
        end
                % wbh draw map with stations and backaz
        subplot(Mphvel,N,length(periods)+1)
        
        ax = worldmap(lalim,lolim);
        set(ax, 'Visible', 'off')
        plotm(pbLat,pbLon,'LineWidth',2,'Color','k') % plate boundaries
        %geoshow(eventcs.stlas,eventcs.stlos,'DisplayType','point','Marker','.','MarkerEdgeColor','b','MarkerSize',12)
        amps = [];
        for ista = 1:length(eventcs.stnms) % loop through stations
            if ismember(ista,badstaids) % don't plot the bad stations
                continue
            end
            amps(ista) = mean(eventcs.autocor(ista).amp);
        end
        scatterm(eventcs.stlas,eventcs.stlos,30,amps,'o','filled')
        caxis([min(amps) max(amps)]);
        h = colorbar;
        ylabel(h,'Amplitude');
        % get event azimuth
        az = azimuth(eventcs.stlas(1),eventcs.stlos(1),evla,evlo);
        [distdeg,az] = distance(eventcs.stlas(1),eventcs.stlos(1),evla,evlo);
        arrlen = 1.4;
        arru = arrlen.*cosd(az);
        arrv = arrlen.*sind(az);
        arrLat = mean(lalim);
        arrLon = mean(lolim);
        quiverm(arrLat,arrLon,arru,arrv,'r')
        azstr = ["az: "+num2str(round(az))+'\circ'];
        textm(lalim(1)+0.5,lolim(1)+0.5,azstr,'FontSize',12) 
        
        MwStr = sprintf('%.2f',eventphv(ip).Mw);
        DistStr = sprintf('%.0f',distdeg);
        
        sgtitle("Phase Velocities for "+eventcs.id+' M'+MwStr+' Dist: '+DistStr+'\circ'); %wbh
        %ofn = [figDir,'PhaseVels_0.25.png'];   % wbh save in event dir
        ofn = [figDirPhv,'/',num2str(round(az)),'_',DistStr,'_M',MwStr,'_PhaseVels_0.25.png'];   % wbh save in BackAz dir

        saveas(gcf,ofn) %wbh
        drawnow;

        % wbh draw ray density
        figure(89)
        clf
        ofn = [figDirPhv,'RayDensity_0.25.png'];
        for ip = 1:length(periods)
            subplot(M,N,ip)
            ax = worldmap(lalim, lolim);
            set(ax, 'Visible', 'off')
            h2 = surfacem(xi,yi,eventphv(ip).raydense);
            title(['Period:', num2str(periods(ip))],'fontsize', 15)
            avgv = nanmean(eventphv(ip).GV(:));
            if isnan(avgv)
                continue;
            end
            maxRayDensity = max(eventphv(ip).raydense,[],'all');
            caxis([0,maxRayDensity])
            colorbar
            colormap summer
            h = colorbar;
            ylabel(h,'Ray Density')
        end
        sgtitle("Ray density for "+eventcs.id);
        saveas(gcf,ofn)

		drawnow;
	end
	matfilename = [eikonl_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat'];
	save(matfilename,'eventphv');
	disp(['Save the result to: ',matfilename])
end % end of loop ie
disp(['Min sta dist of ',num2str(min_stadist_wavelength),' corresponds to:'])
for ip = 1:length(periods)
    disp(['Period: ',num2str(periods(ip)),', dist: ',num2str(ref_phv(ip)*periods(ip)*min_stadist_wavelength),' km.'])
end
