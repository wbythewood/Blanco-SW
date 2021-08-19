% Load all phase velocities and calculate azimuthal anisotropy and isotropic 
% velocity assuming a 1-D structure. 
%
% c(w,theta) = c0(w) * [ 1 + A(w)*cos(2*(theta-phi_A(w))) 
%                          + B(w)*cos(4*(theta-phi_B(w))) ]
%   Russell et al. (2019) JGR DOI:10.1029/2018JB016598
%
% NOTE: For 2-D structure use ./ray_tomo/
%
% https://github.com/jbrussell
clear
setup_parameters;
isoutput_aniso = 1; % write output anisotropy .mat file?
iso_rough = 4.0;
dv_tol = 8;

%======================= PARAMETERS =======================%
InvString = '_Iterate20'; % include leading underscore... 

comp = {parameters.strNAMEcomp};
%xspdir = 'phv_dir';

windir = parameters.winDirName; 
frange = 1./parameters.PeriodRange; 
lalim = parameters.lalim;
lolim = parameters.lolim;

% Quality control parameters:
snr_tol = parameters.tomo_snr_tol; % minimum signal-to-noise
r_tol = parameters.r_tol_min; % [km] minimum separation between stations (should make this number frequency dependent!)
err_tol = parameters.err_tol; % maximum misfit of bessel fit between observed and synthetic
fiterrtol = parameters.fiterrtol; % error allowed in the wavelet fitting
dterrtol = parameters.dterrtol; % largest variance of the inversion error allowed
dep_tol = [0 0]; % [sta1 sta2] OBS depth tolerance
Nwl = parameters.Wavelengths;

% Plotting parameters
ylims_aniso = [-5 5]; %[-3 3];
ylim_p2p = [0 5];

% plate boundaries
mapsDir = [parameters.MapsDir,'PlateBoundaries_NnrMRVL/'];
usgsFN = [parameters.MapsDir,'usgs_plates.txt.gmtdat'];
[pbLat,pbLon] = importPlates(usgsFN);

fastdir = 0; % Expected fast direction for plotting purposes

if comp{1}(1) == 'R'
    ylims = [3.2 4.5];
elseif comp{1}(1) == 'Z'
    ylims = [2.5 4.5];
    %ylims = [2.4 4.5];
elseif comp{1}(1) == 'T'
    ylims = [3.8 4.8];
end
%==========================================================%

%% Load Depths
STAS = stalist;
LATS = stalat;
LONS = stalon;
DEPTHS = staz;

%% Setup Paths

% input path
%XSP_path = [parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];
XSP_path = [parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(Nwl),'wl_phv_dir',InvString,'/'];
% figure output path
%phv_fig_path = [parameters.figpath,windir,'/fullStack/Xsp_anisotropy/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/'];
phv_fig_path = [parameters.XSPfigpath,windir,'/PhV_dir/Ani/snr-',num2str(snr_tol),'minr-',num2str(r_tol),'_fiterr-',num2str(fiterrtol),'_dterr-',num2str(dterrtol),'_err-',num2str(err_tol),'_dvTol-',num2str(dv_tol),InvString,'/'];

if ~exist(phv_fig_path)
    
    mkdir(phv_fig_path);
end

% output path for anisotropy fit
%aniso_path = ['./Xsp/',windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',xspdir,'/azi_aniso_win/'];
aniso_path = [parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(Nwl),'wl_phv_dir',InvString,'/,azi_aniso/'];
if ~exist(aniso_path)   
    mkdir(aniso_path);
end


warning off; %#ok<WNOFF>

stalist = parameters.stalist;
nsta=parameters.nsta; % number of target stations to calculate for

DIRS = dir([XSP_path,'*_xsp.mat']);
nxsp = size(DIRS,1);

%%% --- Loop through XSP files --- %%%
weight_sum = 0;
err_sum = 0;
weight_sum_snr = 0;
snr_sum = 0;
phV = [];
azi = [];
for ixsp=1:nxsp
    
    filename = [XSP_path,DIRS(ixsp).name];

    if ~exist(filename,'file')
        disp(['not exist ',filename])
        continue;
    end

    % LOAD PHV CURVES
    load(filename);

    % SETUP PHV AND AZI ARRAYS
    phV = [phV; xspinfo.r./xspinfo.tw];
    [~,S1az]=distance(xspinfo.lat1,xspinfo.lon1,xspinfo.lat2,xspinfo.lon2);
    if S1az > 180
        S1az = S1az - 360;
    end
    azi = [azi; S1az];
    %if xspinfo.snr >= snr_tol && xspinfo.r >= r_tol
    if xspinfo.snr >= 5 && xspinfo.r >= 100
        sta1 = xspinfo.sta1;
        sta2 = xspinfo.sta2;
        c = xspinfo.r./xspinfo.tw1;
        c_fit = xspinfo.r./xspinfo.tw;
        periods = 1./xspinfo.twloc*2*pi;
        err = xspinfo.sumerr;
        snr = xspinfo.snr;
%        [~,ierr] = min(abs(err-clr_err));
%        [~,isnr] = min(abs(snr-clr_snr));
        [~,S1az]=distance(xspinfo.lat1,xspinfo.lon1,xspinfo.lat2,xspinfo.lon2);

        % Weighted sum error
        weight_sum = weight_sum + c_fit./err;
        err_sum = err_sum + 1/err;

        % Weighted sum snr
        weight_sum_snr = weight_sum_snr + c_fit.*snr;
        snr_sum = snr_sum + snr;

    end
        
end  %end of ixsp
c_weight_avg = weight_sum./err_sum;
c_weight_avg_snr = weight_sum_snr./snr_sum;

phV_QC = [];
azi_QC = [];
err_QC = [];
lats_QC = [];
lons_QC = [];
for ixsp=1:nxsp
    filename = [XSP_path,DIRS(ixsp).name];

    if ~exist(filename,'file')
        disp(['not exist ',filename])
        continue;
    end

    % LOAD PHV CURVES
    load(filename);
    dv_test = (c_fit - iso_rough) ./ iso_rough*100;
    if xspinfo.snr >= snr_tol && xspinfo.r >= r_tol && xspinfo.sumerr <= err_tol ...
            && DEPTHS(strcmp(xspinfo.sta1,STAS)) <= dep_tol(1) ... 
            && DEPTHS(strcmp(xspinfo.sta2,STAS)) <= dep_tol(2) ...
        sta1 = xspinfo.sta1;
        lats_QC = [lats_QC; xspinfo.lat1 xspinfo.lat2];
        lons_QC = [lons_QC; xspinfo.lon1 xspinfo.lon2];
        sta2 = xspinfo.sta2;
        c = xspinfo.r./xspinfo.tw1;
        c_fit = xspinfo.r./xspinfo.tw;
        % wbh add this to remove individual measurements that are problematic
        dv_test = (c_fit - iso_rough) ./ iso_rough*100; % calculate dv (%)
        i_remove = abs(dv_test) <= dv_tol; % find indices for obs that don't meet criterion
        c_fit(i_remove==0) = NaN; %turn them to NaN?
        
        periods = 1./xspinfo.twloc*2*pi;
        err = xspinfo.sumerr;
        snr = xspinfo.snr;
%        [~,ierr] = min(abs(err-clr_err));
%        [~,isnr] = min(abs(snr-clr_snr));
        [~,S1az]=distance(xspinfo.lat1,xspinfo.lon1,xspinfo.lat2,xspinfo.lon2);
        if S1az > 180
            S1az = S1az-360;
        end
        
        % MAKE PHV & AZI ARRAYS
        phV_QC = [phV_QC; c_fit];
        azi_QC = [azi_QC; S1az];
        err_QC = [err_QC; err];
        
        xspinfo.S1az = S1az;
        xspinfo.c_start = c;
        xspinfo.c_fit = c_fit;
        xspinfo.periods = periods;
        aniso.xspinfo(ixsp) = xspinfo;
        
    end
end

% Fit Anisotropy QC
phv_std = std(phV);
for iper = 1:length(periods)
    %varargin{1} = phv_std(iper);
    varargin = sqrt(err_QC);
    %varargin = [];
    %[fitstr(iper), isophv(iper), A_2(iper), phi_2(iper)] = fit_azi_anisotropy(azi,phV(:,iper),varargin,comp{1}(1));
    %[fitstr{iper}, isophv(iper), A2_2(iper), A4_2(iper), phi_2(iper)] = fit_azi_anisotropy2theta4theta(azi_QC,phV_QC(:,iper),comp{1}(1),varargin);
%     [fitstr{iper}, isophv(iper), A2_2(iper), A4_2(iper), phi2_2(iper), phi4_2(iper)] = fit_azi_anisotropy2theta4theta_2(azi_QC,phV_QC(:,iper),comp{1}(1),varargin);
    [fitstr{iper}, isophv(iper), A2_2(iper), A4_2(iper), phi2_2(iper), phi4_2(iper)] = fit_azi_anisotropy2theta(azi_QC,phV_QC(:,iper),comp{1}(1),varargin);
    parastd{iper}=confint(fitstr{iper});
    
    err(iper) = parastd{iper}(2,1) - fitstr{iper}.a;
    err_2p2p(iper) = parastd{iper}(2,2) - fitstr{iper}.d2;
    err_4p2p(iper) = 0;
    err_phi2(iper) = parastd{iper}(2,3) - fitstr{iper}.e2;
    err_phi4(iper) = 0;
    
    % WEIGHTED RMS ERROR
    w = err_QC;
    % 4theta
    dobs_4 = A4_2(iper)*cosd(4*(azi_QC-phi4_2(iper)));
    dpre_4 = (phV_QC(:,iper)-isophv(iper))./isophv(iper) - A2_2(iper)*cosd(2*(azi_QC-phi2_2(iper)));
    wRMS_4A(iper) = sqrt(sum(w.*(dobs_4-dpre_4).^2)/sum(w));
    % 2theta
    dobs_2 = A2_2(iper)*cosd(2*(azi_QC-phi2_2(iper)));
    dpre_2 = (phV_QC(:,iper)-isophv(iper))./isophv(iper) - A4_2(iper)*cosd(4*(azi_QC-phi4_2(iper)));
    wRMS_2A(iper) = sqrt(sum(w.*(dobs_2-dpre_2).^2)/sum(w));
    
end
        
% Phase velocity variations (percent)
c_weight_avg = repmat(c_weight_avg,size(phV_QC,1),1);
%c_perc = (phV_QC-c_weight_avg)./c_weight_avg*100; % from weighted average
isophv_mat = repmat(isophv,size(phV_QC,1),1);
c_perc = (phV_QC-isophv_mat)./isophv_mat*100;

aniso.c_iso = isophv;
aniso.A2 = A2_2;
aniso.A4 = A4_2;
aniso.phi2 = phi2_2;
aniso.phi4 = phi4_2;
aniso.err_c_iso = err;
aniso.err_2A = err_2p2p;
aniso.err_4A = err_4p2p;
aniso.err_phi2 = err_phi2;
aniso.err_phi4 = err_phi4;
aniso.periods = periods;
aniso.fitstr = fitstr;
aniso.wRMS_2A = wRMS_2A;
aniso.wRMS_4A = wRMS_4A;

if isoutput_aniso
    save([aniso_path,'/phv_2theta4theta_wRMS_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'.mat'],'aniso');
end

%% PLOT RAYPATHS
% figure(99); clf;
% cl = lines(5);
% ax = worldmap(lalim, lolim);
% %plot(lons_QC',lats_QC','-k','linewidth',2); hold on;
% %plot(LONS,LATS,'ok','markersize',15,'MarkerFaceColor',cl(2,:),'linewidth',1);
% plotm(lats_QC',lons_QC','Color',[0.3 0.2 0.8],'linewidth',1); hold on;
% plotm(LATS,LONS,'ok','markersize',5,'MarkerFaceColor',[0 0 0],'linewidth',1);
% plotm(pbLat,pbLon,'LineWidth',2,'Color','k') % plate boundaries
% set(gca,'linewidth',2,'fontsize',16,'box','on');
% title('Ray Paths','fontsize',15)
% %grid on;
% xlabel('Longitude'); ylabel('Latitude');

%% Some things for plotting... 
Ipers = 1:length(periods); %[2 3 5 8 10 12];
dimpl = [273   272   761   433];
rowpl = 3;
colpl = 3;
mrksize = 2;
LW = 2;
FS = 15;
symb = 'ok';
mrkclr = [0 0 0];
clr = lines(10);
clr_2theta = clr(2,:);
clr_4theta = clr(1,:);
clr_sum = [0.6 0.6 0.6];
dy_lab = -3.5;

% f3 = figure(3); clf; hold on; set(gcf, 'Color', 'w');
% % set(gcf,'position',[10         248        1203         457]);
% set(gcf,'position',dimpl);
% f5 = figure(5); clf; hold on; set(gcf, 'Color', 'w');
% % set(gcf,'position',[10         248        1203         457]);
% set(gcf,'position',dimpl);


%% PLOT 2-theta (and 4-theta if you have it)
f3 = figure(3);
clf;
ii = 0;
%for iper = Ipers %1:length(periods)
for iper = 1:length(periods)
    ii = ii + 1;
%     subplot(2,4,iper); hold on;
    subplot(rowpl,colpl,ii); hold on;
    a = isophv(iper);
    d2 = A2_2(iper);
    d4 = A4_2(iper);
    e2 = phi2_2(iper);
    e4 = phi4_2(iper);
    if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
        c = 2; % 2 theta
        e_patty = fastdir;
    elseif comp{1}(1) == 'T'
        c = 4; % 4 theta
        e_patty = fastdir-45;
    end
    x = [-180:180];
    % PHV FIT = a*(1+d*cosd(c*(x-e)))
    
%     plot(x,d4*cosd(4*(x-e_patty))*100,'-','color',[.5 .5 .5],'linewidth',LW); %Patty
    % the 4-theta term
    %h1(2) = plot(x,d4*cosd(4*(x-e4))*100,'-','color',clr_4theta,'linewidth',LW);
    % the 2-theta term
    h1(1) = plot(x,d2*cosd(2*(x-e2))*100,'-','color',clr_2theta,'linewidth',LW);
%     h1(5) = plot(x,d2*cosd(2*(x-e2))*100+d4*cosd(4*(x-e4))*100,'-','color',[0 0.7 0.7],'linewidth',LW);
    plot(azi_QC,c_perc(:,iper),symb,'MarkerFaceColor',mrkclr,'markersize',mrksize,'linewidth',1); hold on;
    if iper == 1
%         legend(h1,{'2\theta','4\theta'},'location','northwest');
    end
    if 0 &&  comp{1}(1) == 'T' % plot azimuth bars?
        % fast
        plot([35 35],[-10 10],'--k');
        plot([125 125],[-10 10],'--k');
        plot([215 215]-360,[-10 10],'--k');
        plot([305 305]-360,[-10 10],'--k');

        % slow
        plot([80 80],[-10 10],'-k');
        plot([170 170],[-10 10],'-k');
        plot([260 260]-360,[-10 10],'-k');
        plot([350 350]-360,[-10 10],'-k');
    end
    if comp{1}(1) == 'R' || comp{1}(1) == 'Z'
%         title([num2str(periods(iper)),' s; FAST = ',num2str(e)],'fontsize',FS);
        title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
    elseif comp{1}(1) == 'T'
%         title([num2str(periods(iper)),' s; SLOW = ',num2str(e-45)],'fontsize',FS);
        title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
    end
    if ii == 9
        xlabel('Azimuth (degrees)','fontsize',FS);
    end
    if ii == 1
        yl = ylabel('\delta{c}/c (%)','fontsize',FS);
        yl.Position = [yl.Position(1) yl.Position(2)+dy_lab yl.Position(3:end)];
    end
    set(gca,'fontsize',FS,'linewidth',1.5);
    xlim([-180 180]);
    %ylim(ylims_aniso);
    %ylim([3.8 4.8]);
    %box on;
    
    err(iper) = parastd{iper}(2,1) - fitstr{iper}.a;
end
sgtitle('2-Theta variations') 

%% Plot 2-theta + 4-theta
% f5 = figure(5);
% ii = 0;
% for iper = Ipers
%     ii = ii + 1;
% %     subplot(2,4,iper); hold on;
%     subplot(rowpl,colpl,ii); hold on;
%     a = isophv(iper);
%     d2 = A2_2(iper);
%     d4 = A4_2(iper);
%     e2 = phi2_2(iper);
%     e4 = phi4_2(iper);
%     if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
%         c = 2; % 2 theta
%         e_patty = fastdir;
%     elseif comp{1}(1) == 'T'
%         c = 4; % 4 theta
%         e_patty = fastdir-45;
%     end
%     x = [-180:180];
%     % PHV FIT = a*(1+d*cosd(c*(x-e)))
%     
%     h2(3) = plot(x,d2*cosd(2*(x-e2))*100+d4*cosd(4*(x-e4))*100,'-','color',clr_sum,'linewidth',LW);
%     plot(azi_QC,c_perc(:,iper),symb,'MarkerFaceColor',mrkclr,'markersize',mrksize,'linewidth',1); hold on;
%     if iper == 1
% %         legend(h2,{'2\theta + 4\theta'},'location','northwest');
%     end
%     if comp{1}(1) == 'R' || comp{1}(1) == 'Z'
% %         title([num2str(periods(iper)),' s; FAST = ',num2str(e)],'fontsize',FS);
%         title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
%     elseif comp{1}(1) == 'T'
% %         title([num2str(periods(iper)),' s; SLOW = ',num2str(e-45)],'fontsize',FS);
%         title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
%     end
%     if ii == 5
%         xlabel('Azimuth (degrees)','fontsize',FS);
%     end
%     if ii == 1
%         yl = ylabel('\delta{c}/c (%)','fontsize',FS);
%         yl.Position = [yl.Position(1) yl.Position(2)+dy_lab yl.Position(3:end)];
%     end
%     set(gca,'fontsize',FS,'linewidth',1.5);
%     xlim([-180 180]);
%     %ylim(ylims_aniso);
%     %ylim([3.8 4.8]);
%     %box on;
%     
%     err(iper) = parastd{iper}(2,1) - fitstr{iper}.a;
% end
% sgtitle('2-Theta + 4-Theta variations') 

%% Plot Phase Velocities

f54 = figure(54); 
clf;
h54(1) = plot(periods,c_weight_avg(1,:),'-ok','linewidth',2); hold on;
%plot(periods,isophv,'-or','linewidth',2);
h54(2) = errorbar(periods,isophv,err,'-or','linewidth',2);
h54(3) = plot(periods,xspinfo.c_start,'-o','linewidth',2,'color',[0.5 0.5 0.5]);

ylim(ylims);
xlim([1/frange(1) 1/frange(2)]);
xlabel('Period (s)','fontsize',15);
ylabel('Phase Velocity (km/s)','fontsize',15);
set(gca,'fontsize',12);
legend(h54,{'PHV_{avg}','PHV_{iso}','PHV_{start}'},'fontsize',12,'location','southeast');

%%
% PLOT Amplitude and Fast direction
f4 = figure(4); clf;
% set(gcf,'position',[4         325        1239         380],'color','w');
set(gcf,'position',[6   220   490   485]);
for iper = 1:length(periods)
%     err_p2p(iper) = parastd{iper}(2,2) - fitstr{iper}.d;
    err_2p2p(iper) = parastd{iper}(2,2) - fitstr{iper}.d2;
    err_4p2p(iper) = 0;
    err_phi2(iper) = parastd{iper}(2,3) - fitstr{iper}.e2;
    err_phi4(iper) = 0;
end
% peak-to-peak
subplot(2,1,1); hold on;
% plot(periods,ones(size(periods))*4.2,'--k','linewidth',2);
if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
    h3(1) = errorbar(periods,A2_2*2*100,wRMS_2A*100,'-o','color',clr_2theta,'linewidth',2);
else
    h3(2) = errorbar(periods,A4_2*2*100,wRMS_4A*100,'-o','color',clr_4theta,'linewidth',2);
    h3(1) = errorbar(periods,A2_2*2*100,wRMS_2A*100,'-o','color',clr_2theta,'linewidth',2);
end
% xlim([3.5 10.5]);
%xlim([1/frange(2) 1/frange(1)]);
xlim([1/frange(1) 1/frange(2)]);
%ylim(ylim_p2p);
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',FS);
xlabel('Period (s)','fontsize',FS);
ylabel('Peak-to-peak amp (%)','fontsize',FS);
% legend(h3,{'2\theta','4\theta'},'location','northwest');

% Azimuth
subplot(2,1,2); hold on;
if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
    for iper = 1:length(periods)
        phi2_vec(1) = phi2_2(iper);
        phi2_vec(2) = phi2_2(iper)+180;
        phi2_vec(3) = phi2_2(iper)-180;
        [~, I] = min(abs(phi2_vec-fastdir));
        phi2_2(iper) = phi2_vec(I);


        phi4_vec(1) = phi4_2(iper);
        phi4_vec(2) = phi4_2(iper)+90;
        phi4_vec(3) = phi4_2(iper)+180;
        phi4_vec(4) = phi4_2(iper)+270;
        phi4_vec(5) = phi4_2(iper)-90;
        phi4_vec(6) = phi4_2(iper)-180;
        phi4_vec(7) = phi4_2(iper)-270;
        [~, I] = min(abs(phi4_vec-fastdir));
        phi4_2(iper) = phi4_vec(I);
    end
    % for horizontal lines...
    %plot(periods,ones(size(periods))*fastdir,'--k','linewidth',2);
    %plot(periods,ones(size(periods))*(fastdir+45),'--k','linewidth',2);
    %plot(periods,ones(size(periods))*(fastdir+90),'--k','linewidth',2);
    %
%     errorbar(periods,phi4_2,err_phi4,'-ob','linewidth',2);
    errorbar(periods,phi2_2,err_phi2,'-o','color',clr_2theta,'linewidth',2);
    ylabel('Fast Direction (%)','fontsize',FS);
    %ylim([0 180]);
end
if comp{1}(1) == 'T'
    for iper = 1:length(periods)
        phi2_vec(1) = phi2_2(iper);
        phi2_vec(2) = phi2_2(iper)+180;
        phi2_vec(3) = phi2_2(iper)-180;
        [~, I] = min(abs(phi2_vec-fastdir+90));
        phi2_2(iper) = phi2_vec(I);
        if phi2_2(iper) < fastdir
            phi2_2(iper) = phi2_2(iper)+180;
        end


        phi4_vec(1) = phi4_2(iper);
        phi4_vec(2) = phi4_2(iper)+90;
        phi4_vec(3) = phi4_2(iper)+180;
        phi4_vec(4) = phi4_2(iper)+270;
        phi4_vec(5) = phi4_2(iper)-90;
        phi4_vec(6) = phi4_2(iper)-180;
        phi4_vec(7) = phi4_2(iper)-270;
        [~, I] = min(abs(phi4_vec-fastdir+45));
        phi4_2(iper) = phi4_vec(I);
        if phi4_2(iper) < fastdir
            phi4_2(iper) = phi4_2(iper)+90;
        end
    end
    
    plot(periods,ones(size(periods))*fastdir,'--k','linewidth',2);
    plot(periods,ones(size(periods))*(fastdir+45),'--k','linewidth',2);
    plot(periods,ones(size(periods))*(fastdir+90),'--k','linewidth',2);
%     plot(periods,ones(size(periods))*(78+45),'--k','linewidth',2);
%     plot(periods,ones(size(periods))*(78+90),'--k','linewidth',2);

%     errorbar(periods,phi4_2+45,err_phi4,'-ob','linewidth',2);
%     errorbar(periods,phi2_2+90,err_phi2,'-o','color',clr_2theta,'linewidth',2);
    errorbar(periods,phi4_2,err_phi4,'-o','color',clr_4theta,'linewidth',2);
    errorbar(periods,phi2_2,err_phi2,'-o','color',clr_2theta,'linewidth',2);
    ylabel('Fast Direction (\circ)','fontsize',FS);
    ylim([50 200]);
end
% xlim([3.5 10.5]);
xlim([1/frange(1) 1/frange(2)]);
% xlim([4.5 10.1]);
set(gca,'linewidth',1.5,'xminortick','on','yminortick','on','fontsize',FS);
xlabel('Period (s)','fontsize',FS,'linewidth',1.5);

%% Subtract 4 theta signal
% f6 = figure(6); clf; hold on; set(gcf, 'Color', 'w');
% set(gcf,'position',dimpl);
% ii = 0;
% for iper = Ipers
%     ii = ii + 1;
% %     subplot(2,4,iper); hold on;
%     subplot(rowpl,colpl,ii); hold on;
%     a = isophv(iper);
%     d2 = A2_2(iper);
%     d4 = A4_2(iper);
%     e2 = phi2_2(iper);
%     e4 = phi4_2(iper);
%     if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
%         c = 2; % 2 theta
%         e_patty = fastdir;
%     elseif comp{1}(1) == 'T'
%         c = 4; % 4 theta
%         e_patty = fastdir-45;
%     end
%     x = [-180:180];
%     % PHV FIT = a*(1+d*cosd(c*(x-e)))
%     
% %     patch([x fliplr(x)],[(d2+wRMS_2A(iper))*cosd(2*(x-e2))*100 fliplr((d2-wRMS_2A(iper))*cosd(2*(x-e2))*100)],[0.8 0.8 0.8],'linestyle','none');
%     h1(1) = plot(x,d2*cosd(2*(x-e2))*100,'-','color',clr_2theta,'linewidth',LW);
%     plot(azi_QC,c_perc(:,iper)-d4*cosd(4*(azi_QC-e4))*100,symb,'MarkerFaceColor',mrkclr,'markersize',mrksize,'linewidth',1); hold on;
% 
%     
%     if iper == 1
% %         legend(h1,{'2\theta','4\theta'},'location','northwest');
%     end
%     if comp{1}(1) == 'R' || comp{1}(1) == 'Z'
% %         title([num2str(periods(iper)),' s; FAST = ',num2str(e)],'fontsize',FS);
%         title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
%     elseif comp{1}(1) == 'T'
% %         title([num2str(periods(iper)),' s; SLOW = ',num2str(e-45)],'fontsize',FS);
%         title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
%     end
%     if ii == 5
%         xlabel('Azimuth (degrees)','fontsize',FS);
%     end
%     if ii == 1
%         yl = ylabel('\delta{c}/c (%)','fontsize',FS);
%         yl.Position = [yl.Position(1) yl.Position(2)+dy_lab yl.Position(3:end)];
%     end
%     set(gca,'fontsize',FS,'LineWidth',1.5);
%     xlim([-180 180]);
%     %ylim(ylims_aniso);
%     %ylim([3.8 4.8]);
%     %box on;
%     
%     err(iper) = parastd{iper}(2,1) - fitstr{iper}.a;
% end

%% Subtract 2 theta signal
f7 = figure(7); clf; hold on; set(gcf, 'Color', 'w');
set(gcf,'position',dimpl);
ii = 0;
for iper = Ipers
    ii = ii + 1;
%     subplot(2,4,iper); hold on;
    subplot(rowpl,colpl,ii); hold on;
    a = isophv(iper);
    d2 = A2_2(iper);
    d4 = A4_2(iper);
    e2 = phi2_2(iper);
    e4 = phi4_2(iper);
    if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
        c = 2; % 2 theta
        e_patty = fastdir;
    elseif comp{1}(1) == 'T'
        c = 4; % 4 theta
        e_patty = fastdir-45;
    end
    x = [-180:180];
    % PHV FIT = a*(1+d*cosd(c*(x-e)))
    
%     patch([x fliplr(x)],[(d4+wRMS_4A(iper))*cosd(4*(x-e4))*100 fliplr((d4-wRMS_4A(iper))*cosd(4*(x-e4))*100)],[0.8 0.8 0.8],'linestyle','none');
    h1(2) = plot(x,d4*cosd(4*(x-e4))*100,'-','color',clr_4theta,'linewidth',LW);
    plot(azi_QC,c_perc(:,iper)-d2*cosd(2*(azi_QC-e2))*100,symb,'MarkerFaceColor',mrkclr,'markersize',mrksize,'linewidth',1); hold on;
    if iper == 1
%         legend(h1,{'2\theta','4\theta'},'location','northwest');
    end
    if comp{1}(1) == 'R' || comp{1}(1) == 'Z'
%         title([num2str(periods(iper)),' s; FAST = ',num2str(e)],'fontsize',FS);
        title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
    elseif comp{1}(1) == 'T'
%         title([num2str(periods(iper)),' s; SLOW = ',num2str(e-45)],'fontsize',FS);
        title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
    end
    if ii == 5
        xlabel('Azimuth (degrees)','fontsize',FS);
    end
    if ii == 1
        yl = ylabel('\delta{c}/c (%)','fontsize',FS);
        yl.Position = [yl.Position(1) yl.Position(2)+dy_lab yl.Position(3:end)];
    end
    set(gca,'fontsize',FS,'linewidth',1.5);
    xlim([-180 180]);
    %ylim(ylims_aniso);
    %ylim([3.8 4.8]);
    %box on;
    
end
sgtitle('Subtracting the 2-Theta trend') 

%% Subtract 2theta and 4theta signal
% f8 = figure(8); clf; hold on; set(gcf, 'Color', 'w');
% set(gcf,'position',dimpl);
% ii = 0;
% for iper = Ipers
%     ii = ii + 1;
% %     subplot(2,4,iper); hold on;
%     subplot(rowpl,colpl,ii); hold on;
%     a = isophv(iper);
%     d2 = A2_2(iper);
%     d4 = A4_2(iper);
%     e2 = phi2_2(iper);
%     e4 = phi4_2(iper);
%     if comp{1}(1) == 'Z' || comp{1}(1) == 'R'
%         c = 2; % 2 theta
%         e_patty = fastdir;
%     elseif comp{1}(1) == 'T'
%         c = 4; % 4 theta
%         e_patty = fastdir-45;
%     end
%     x = [-180:180];
%     % PHV FIT = a*(1+d*cosd(c*(x-e)))
%     
%     residual = c_perc(:,iper)-d2*cosd(2*(azi_QC-e2))*100-d4*cosd(4*(azi_QC-e4))*100;
%     plot(azi_QC,residual,symb,'MarkerFaceColor',mrkclr,'markersize',mrksize,'linewidth',1); hold on;
%     plot([-180 180],[1 1]*rms(residual),'--k');
%     plot([-180 180],[-1 -1]*rms(residual),'--k');
%     if comp{1}(1) == 'R' || comp{1}(1) == 'Z'
%         title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
%     elseif comp{1}(1) == 'T'
%         title([num2str(periods(iper),'%0.1f'),' s'],'fontsize',30);
%     end
%     if ii == 5
%         xlabel('Azimuth (degrees)','fontsize',FS);
%     end
%     if ii == 1
%         yl = ylabel('\delta{c}/c (%)','fontsize',FS);
%         yl.Position = [yl.Position(1) yl.Position(2)+dy_lab yl.Position(3:end)];
%     end
%     set(gca,'fontsize',FS,'linewidth',1.5);
%     xlim([-180 180]);
%     ylim(ylims_aniso);
% 
% end

%%
psfile1 = [phv_fig_path,'xsp_fitanisotropy2theta4theta_2_wRMS_win_',comp{1}(1),'_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'_Ani-data.pdf'];
%print('-dpsc2',psfile);
%save2pdf(psfile1,f3,1000);
saveas(f3,psfile1);
psfile2 = [phv_fig_path,'xsp_fitanisotropy2theta4theta_2_wRMS_win_',comp{1}(1),'_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'_A-Phi.pdf'];
%print('-dpsc2',psfile);
%save2pdf(psfile2,f4,1000);
saveas(f4,psfile2);
psfile3 = [phv_fig_path,'xsp_fitanisotropy2theta4theta_2_wRMS_win_',comp{1}(1),'_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'_sum_ES17.pdf'];
%print('-dpsc2',psfile);
%save2pdf(psfile3,f5,1000);
%saveas(f5,psfile3);

psfile6 = [phv_fig_path,'xsp_fitanisotropy2theta4theta_2_wRMS_win_',comp{1}(1),'_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'_2theta_ES17.pdf'];
%print('-dpsc2',psfile);
%save2pdf(psfile6,f6,1000);
%saveas(f6,psfile6);

psfile7 = [phv_fig_path,'xsp_fitanisotropy2theta4theta_2_wRMS_win_',comp{1}(1),'_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'_4theta_ES17.pdf'];
%print('-dpsc2',psfile);
%save2pdf(psfile7,f7,1000);
saveas(f7,psfile7);

psfile8 = [phv_fig_path,'xsp_fitanisotropy2theta4theta_2_wRMS_win_',comp{1}(1),'_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'_Resid.pdf'];
%print('-dpsc2',psfile);
%save2pdf(psfile8,f8,1000);
%saveas(f8,psfile8);

psfile54 = [phv_fig_path,'xsp_fitanisotropy2theta4theta_2_wRMS_win_',comp{1}(1),'_SNRtol',num2str(snr_tol),'_disttol',num2str(r_tol),'_errtol',num2str(err_tol),'_PhVel.pdf'];
saveas(f54,psfile54);
