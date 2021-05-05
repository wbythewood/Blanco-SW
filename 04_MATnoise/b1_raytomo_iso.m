% Script to do the ray theory tomography based on the ambient noise measurement. 
% Solves for azimuthal anisotropy and isotropic velocity on a 2D grid.
%
% phv(theta,freq) = phv_iso(freq) + Ac2(freq)*cos(2*theta) + As2(freq)*sin(2*theta)
%                                 + Ac4(freq)*cos(4*theta) + As4(freq)*sin(4*theta)
%
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012
% JBR: Modified in 2018 to include azimuthal anisotropy
% https://github.com/jbrussell

% wbh modified to include only isotropic velocity
clear; close all;

%%
% Save results?
isoutput = 1;
savefile = ['test'];

%Zero Crossing Options
% 0 is only bessel, 1 is ZC as first pass, then bessel fit
isZC = 1;
%======================= PARAMETERS =======================%
setup_parameters;
comp = {parameters.strNAMEcomp};
windir = parameters.winDirName; 
frange = 1./parameters.PeriodRange; 
Nwl = parameters.Wavelengths;

average_vel = 3.8; % [km/s] For calculating wavelength for determining r_tol_min

% QC parameters
snr_tol = parameters.tomo_snr_tol;
is_rtolmin_wavelength = parameters.is_rtolmin_wavelength; 
wl_fac = parameters.wl_fac;
r_tol_min = parameters.r_tol_min;
r_tol_max = parameters.r_tol_max;
err_tol = parameters.err_tol; 
fastdir = 78; % Fast direction for azimuthal anisotropy (only for plotting purposes);
iscompare_aniso = 0; % compare to old anisotropic measurements

%==========================================================%
%%
% Load color scale
load ../seiscmap.mat
load ../roma.mat
% Load anisotropy data (from old inversion)
if iscompare_aniso
    load(['./aniso_DATA/phv_dir/',aniso_data]);
end

% Set up geometry parameters
setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
station_list = parameters.StaListFile;

% Load station info
[sta.nw sta.name, sta.lat, sta.lon, sta.dep] = textread(station_list,'%s %s %f %f %f');

fiterrtol = parameters.fiterrtol;
maxerrweight = parameters.maxerrweight;
polyfit_dt_err = parameters.polyfit_dt_err;
smweight0 = parameters.smweight0;
dterrtol = parameters.dterrtol;
raydensetol = parameters.raydensetol;
r = parameters.r;
%phv_fig_path = [parameters.figpath,windir,'/,PhV_dir/',num2str(1/frange(1)),'_',num2str(1/frange(2)),'s/raytomo_azi2theta_2D/'];
if isZC == 0
    phv_fig_path = [parameters.XSPfigpath,windir,'/PhV_dir/Iso/sm-',num2str(smweight0),'_grid-',num2str(gridsize),'_fiterr-',num2str(fiterrtol),'_dterr-',num2str(dterrtol),'/'];
elseif isZC == 1
    phv_fig_path = [parameters.XSPfigpath,windir,'/PhV_dir/Iso/sm-',num2str(smweight0),'_grid-',num2str(gridsize),'_fiterr-',num2str(fiterrtol),'_dterr-',num2str(dterrtol),'_ZC-Bessel/'];
end

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx = length(xnode);
Ny = length(ynode);

% figure output path
%phv_fig_path = ['./figs/',windir,'/fullStack/raytomo_azi2theta_2D/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/',stafile,'/'];
%phv_fig_path = ['./figs/',windir,'/fullStack/raytomo_azi2theta_2D/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_phv_dir/',station_list,'/'];
if ~exist(phv_fig_path)    
    mkdir(phv_fig_path);
end

% read in bad station list, if existed
if exist('badsta.lst')
    badstnms = textread('badsta.lst','%s');
    badstaids = find(ismember({stainfo.staname},badstnms));
    disp('Found Bad stations:')
    disp(badstnms)
end

% plate boundaries
mapsDir = [parameters.MapsDir,'PlateBoundaries_NnrMRVL/'];
usgsFN = [parameters.MapsDir,'usgs_plates.txt.gmtdat'];
[pbLat,pbLon] = importPlates(usgsFN);

%% Set up initial smoothing kernels (second derivative)
% Isotropic smoothing kernels
F_iso = smooth_kernel_build(xnode, ynode, Nx*Ny);

F = F_iso;

%%
% Initialize the xsp structure
if isZC == 0
    Xsp_path = [parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(Nwl),'wl_phv_dir/'];
elseif isZC == 1
    Xsp_path = [parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(Nwl),'wl_phv_dir_ZC-Bessel/'];
end
xspfiles = dir([Xsp_path,'*_xsp.mat']);

disp('Looking at Xsp Files')
for ixsp = 1:length(xspfiles)
    
    temp = load([Xsp_path,xspfiles(ixsp).name]);
    xspinfo = temp.xspinfo;
    %wbh why is waxis only retreived when ixsp is 1? if it changes from
    %file to file, there is an error later... try to retrieve it each time
    waxis = temp.waxis;
    
    if ixsp ==1
        Tperiods = (2*pi)./temp.twloc;
        %waxis = temp.waxis;
        twloc = temp.twloc;
        xspinfo.isgood = zeros(size(Tperiods));
        xspsum = xspinfo;
        wavelength = average_vel*Tperiods;
    else
        xspinfo.isgood = zeros(size(Tperiods));
        xspsum = [xspsum;xspinfo];
    end
    clear temp

    
    % 	xspinfo(ixsp).isgood = 0;
%     if xspsum(ixsp).sumerr < errlevel ...
%             && xspsum(ixsp).snr > snrtol && xspsum(ixsp).coherenum > mincoherenum
%         xspsum(ixsp).isgood = 1;
%     end

    for ip = 1:length(Tperiods)
        if ~is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= r_tol_min && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        elseif  is_rtolmin_wavelength && xspinfo.snr >= snr_tol && xspinfo.r >= wavelength(ip)*wl_fac && xspinfo.r <= r_tol_max && xspinfo.sumerr <= err_tol
            xspsum(ixsp).isgood(ip) = 1;
        end

        if isempty(cell2mat(strfind(sta.name,xspsum(ixsp).sta1))) || isempty(cell2mat(strfind(sta.name,xspsum(ixsp).sta2)))
            xspsum(ixsp).isgood(ip) = 0;
        end
    end
    
    if rem(ixsp,500)==0
        disp(['Looking at #',num2str(ixsp),' of ',num2str(length(xspfiles))])
    end
end % end of loop ixsp'


% Loop through periods
for ip=1:length(Tperiods)
    disp(' ');
    disp(['Inversing Period: ',num2str(Tperiods(ip))]);
    clear rays dt fiterr mat phaseg err raydense dist azi mat_azi phv
    raynum = 0;

    for ixsp = 1:length(xspsum)
        if xspsum(ixsp).isgood(ip) ==0
            continue;
        end
%         if xspsum(ixsp).r > refv*Tperiods(ip)*distrange(2)...
%                 || xspsum(ixsp).r < refv*Tperiods(ip)*distrange(1)
%             continue;
%         end
        
        raynum = raynum+1;
        rays(raynum,1) = xspsum(ixsp).lat1;
        rays(raynum,2) = xspsum(ixsp).lon1;
        rays(raynum,3) = xspsum(ixsp).lat2;
        rays(raynum,4) = xspsum(ixsp).lon2;
        
        % dist(raynum) = deg2km(distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4)));
        dist(raynum) = distance(rays(raynum,1),rays(raynum,2),rays(raynum,3),rays(raynum,4),referenceEllipsoid('GRS80'))/1000;
        if isZC == 0
            dt(raynum) = xspsum(ixsp).tw(ip);
        elseif isZC == 1
            dt(raynum) = xspsum(ixsp).tw_bfit(ip);
        end
        phv(raynum) = dist(raynum)./dt(raynum);
        
        dep1 = sta.dep(strcmp(xspsum(raynum).sta1,sta.name));
        dep2 = sta.dep(strcmp(xspsum(raynum).sta2,sta.name));
        dep(raynum) = mean([dep1 dep2]);
        
        err = smooth((abs(xspsum(ixsp).err)./mean(abs(xspsum(ixsp).xsp))).^2,round(length(waxis)/length(twloc)));
        fiterr(raynum) = interp1(waxis(:),err(:),twloc(ip)); 
        % Fix the fact that last period always breaks (JBR 9/29/17)
        if isnan(fiterr(raynum))
            [~,I] = min(abs(twloc(ip)-waxis(:)));
            fiterr(raynum) = err(I);
        end
        csnum(raynum) = xspsum(ixsp).coherenum;
        snr(raynum) = xspsum(ixsp).snr;
        errays(raynum,1) = xspsum(ixsp).lat1;
        errays(raynum,2) = xspsum(ixsp).lon1;
        errays(raynum,3) = xspsum(ixsp).lat2;
        errays(raynum,4) = xspsum(ixsp).lon2; 
        errays(raynum,5) = fiterr(raynum);
        
    end
    if size(dt,1) ~=raynum
        dt = dt';
    end
    
    % Building the isotropic data kernel
    disp('Start building the kernel');
    tic
    mat_iso=ray_kernel_build(rays,xnode,ynode);   
    toc
    
    % JBR - Combine isotropic and anisotropic
    mat = mat_iso;
    
    % Calculate the weighting matrix
    W = sparse(length(dt),length(dt));
    for i=1:length(dt)
        W(i,i)=1./fiterr(i);
    end
    ind = find(W > maxerrweight);
    W(ind) = maxerrweight;
    ind = find(W < 1/fiterrtol);
    W(ind) = 0;
    for i=1:length(dt)
        W(i,i)=W(i,i).*(csnum(i).^0.5);
    end
    para = polyfit(dist(:),dt,1);
    polyerr = polyval(para,dist(:)) - dt;
    errind = find(abs(polyerr) > polyfit_dt_err);
    for i = errind
        W(i,i) = 0;
    end
    
    % calculate the smoothing weight
    NR=norm(F,1);
    NA=norm(W*mat,1);
    smweight = smweight0*NA/NR;

    disp('start inversion...');

    A=[W*mat; smweight*F];
    rhs=[W*dt; zeros(size(F,1),1)];

    phaseg=(A'*A)\(A'*rhs);   
    
    % Iteratively down weight the measurement with high error
    niter=1;
    
    while niter < 2
        niter=niter+1;
        err = mat*phaseg - dt;

        stderr=std(err);
        if stderr > dterrtol
            stderr = dterrtol;
        end
        ind = find(diag(W)==0);
        disp('Before iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        for i=1:length(err)
            if abs(err(i)) > 2*stderr
                W(i,i)=0;
            end
        end
        ind = find(diag(W)==0);
        disp('After iter:');
        disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
        disp(['Bad Measurement Number: ', num2str(length(ind))]);
        
        % Rescale the smooth kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;
        
        % Invert

        A=[W*mat; smweight*F];
        rhs=[W*dt; zeros(size(F,1),1)];

        phaseg=(A'*A)\(A'*rhs);

        
    end
    
    % Isotropic phase velocity
    phv_iso = dist'./(mat_iso*phaseg(1:Nx*Ny));
        
    Igood = find(diag(W)~=0);
    mat_good = mat(Igood,:);
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            %raydense(i,j) = sum(mat(:,n));
            raydense(i,j) = sum(mat_good(:,n));
            if raydense(i,j) < raydensetol
                phaseg(n)=NaN;
            end
        end
    end
    
    % Convert into phase velocity
    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            GV(i,j)= 1./phaseg(n);
        end
    end
    
    % JBR - Get Azimuthal coefficients from phaseg (s/km)
    phv_av = nanmean(GV(:));
%     phv_av = nanmedian(GV(:));
    phv_avstd = nanstd(GV(:));
    slow_av = 1/phv_av;
%     slow_av = 1/mean(phv_iso);

    raytomo(ip).GV = GV;
    raytomo(ip).mat = mat;
    raytomo(ip).raydense = raydense;
    raytomo(ip).period = Tperiods(ip);
    raytomo(ip).w = diag(W);
    raytomo(ip).err = err;
    raytomo(ip).rays = rays;
    raytomo(ip).fiterr = fiterr;
    raytomo(ip).dt = dt;
    raytomo(ip).smweight0 = smweight0;
    %JBR    
    raytomo(ip).phv_iso = phv_iso;    
    raytomo(ip).phv_av = phv_av;
    raytomo(ip).phv_avstd = phv_avstd;
    raytomo(ip).phv = phv;

    
    
    if 0
        f1 = figure(1);
        clf
        ax = worldmap(lalim, lolim);
        set(ax, 'Visible', 'off')
        surfacem(xi,yi,raytomo(ip).GV);
%         drawlocal
        title([num2str(Tperiods(ip))],'fontsize',15)
        avgv = nanmean(raytomo(ip).GV(:));
        caxis([avgv*(1-r) avgv*(1+r)])
        colorbar
        colormap(seiscmap)
        ofn = [phv_fig_path,'.pdf'];
        
%         pause;
    end
    
end % end of period loop

lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
% isoutput = 1;
if isoutput
    save(savefile,'raytomo','xnode','ynode');
    save('coor.mat','xi','yi','xnode','ynode','gridsize','lalim','lolim');
end
for iper = 1:length(Tperiods)
    phv_av_rt(iper) = raytomo(iper).phv_av;
    phv_avstd_rt(iper) = raytomo(iper).phv_avstd;
end
%% Phase Velocity Maps (km/s)

Mp = 3; Np = 3;
fig17 = figure(17);
set(gcf,'position',[1    1   1244   704]);
clf
for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    set(gcf,'color',[0.9 0.9 0.9])
    avgv = nanmean(raytomo(ip).GV(:));
    levels = linspace(avgv*(1-r), avgv*(1+r),10);
    contourfm(xi,yi,raytomo(ip).GV,levels);
    plotm(pbLat,pbLon,'LineWidth',2,'Color','k') % plate boundaries
    title([num2str(Tperiods(ip))],'fontsize',15)
    caxis([avgv*(1-r) avgv*(1+r)])
    cb = colorbar;
    ylabel(cb,'Vs (km/s)');
    colormap(roma)
    
    hold on;
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
end
ofn = [phv_fig_path,'phv_km-s.png'];
saveas(fig17,ofn);

%% Phase Velocity Maps (%)

fig19 = figure(19);
set(gcf,'position',[1    1   1244   704]);
clf
vperc = [-r r];
for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    set(gcf,'color',[0.9 0.9 0.9])
    avgv = nanmean(raytomo(ip).GV(:));
    resid = (raytomo(ip).GV-avgv)./avgv;
    levels = linspace(vperc(1),vperc(2),10)*100;
    contourfm(xi,yi,resid*100,levels);
    plotm(pbLat,pbLon,'LineWidth',2,'Color','k') % plate boundaries
    title([num2str(Tperiods(ip))],'fontsize',15)
    caxis([min(levels) max(levels)])
    cb = colorbar;
    ylabel(cb,'dVs (%)');
    colormap(roma)
    
    hold on;
    plotm(sta.lat,sta.lon,'ok','markerfacecolor',[0 0 0]);
end
ofn = [phv_fig_path,'phv_percent.png'];
saveas(fig19,ofn);

%% RAY DENSITY
fig18 = figure(18);
set(gcf,'position',[1    1   1244   704]);
clf

for ip=1:length(Tperiods)
subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).raydense);
%     drawlocal
    plotm(pbLat,pbLon,'LineWidth',2,'Color','k') % plate boundaries
    title([num2str(Tperiods(ip))],'fontsize',15)
    cb = colorbar;
    ylabel(cb,'Ray density');
    colormap(flip(hot));
    caxis([0 500])
    caxis([0 max(raytomo(ip).raydense(:))])
end
ofn = [phv_fig_path,'ray_density.png'];
saveas(fig18,ofn);

%% ERRORS ON XSP
fig16 = figure(16);
set(gcf,'position',[1    1   1244   704]);
clf

for ip=1:length(Tperiods)
    subplot(Mp,Np,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    clear rays
    rays = raytomo(ip).rays;
%    surfacem(xi,yi,raytomo(ip).err);
%    drawpng
    scatterm((rays(:,1)+rays(:,3))./2,(rays(:,2)+rays(:,4))./2,30,raytomo(ip).fiterr,'filled')
%    drawlocal
    plotm(pbLat,pbLon,'LineWidth',2,'Color','k') % plate boundaries
    title([num2str(Tperiods(ip))],'fontsize',15)
    colorbar
    caxis([ 0 5])

end
ofn = [phv_fig_path,'XSP-err.png'];
saveas(fig16,ofn);

%% Plot phase velocities
for ip = 1:length(Tperiods)
    avgv(ip) = nanmean(raytomo(ip).GV(:));
    avgv_std(ip) = nanstd(raytomo(ip).GV(:));
end

fig2 = figure(2);clf;
set(gcf,'position',[227 385 967 324]);

hold on; box on;
try 
    plot(Tperiods,xspinfo.c_start,'ok','linewidth',2);
catch
    display('No starting data')
end
errorbar(Tperiods,phv_av_rt,phv_avstd_rt*2,'-r','linewidth',2);
% errorbar(Tperiods,mean([vertcat(raytomo(:).phv)],2),std([vertcat(raytomo(:).phv)],0,2)*2,'-b','linewidth',2);
title('Isotropic Phase Velocity');
xlabel('Period (s)','fontsize',16);
ylabel('Phase Velocity (km/s)','fontsize',16);
set(gca,'fontsize',16,'linewidth',1.5);
legend({'Starting','Raytomo Avg.'},'location','southeast','fontsize',12,'box','off');
xlimvals = [1./frange]; %wbh make sure xlimvals is increasing
if xlimvals(1) > xlimvals(2)
    xlimvals = flip(xlimvals);
end
xlim(xlimvals);
if comp{1}(1) == 'Z'
    %ylim([3.4 4.3]);
    ylim([3.0 5.0]);
elseif comp{1}(1) == 'T'
    ylim([3.8 4.7]);
end

ofn = [phv_fig_path,'IsoPhaseVelocity.png'];
saveas(fig2,ofn);
%% Save parameters

%also save all the parameters in the fig file so it can easily be figured
%out later...
ofn = [phv_fig_path,'parameters.mat'];
save(ofn,'parameters')