% Extract phase velocity dispersion between station pairs by fitting J0 bessel 
% function to real(ccf)
% Uses cross spectral fitting technique of Menke & Jin (2015) BSSA 
% DOI:10.1785/0120140245
%
% Define own starting phase velocity dispersion c manually or using 
% functions/calc_Rayleigh_disp for a simple layered model (does not work for 
% models with a water column).
%
% https://github.com/jbrussell
clear
close all;

global tN
global waxis
global twloc
global weight
setup_parameters;

%======================= PARAMETERS =======================%
is_resume = 1; % Resume from last processed file (1) or overwrite (0)
isoutput = 1; % Save *.mat file with results?
npts_smooth = 1; % 1 = no smoothing

comp = {parameters.strNAMEcomp};
windir = parameters.winDirName; 
figDir = parameters.XSPfigpath;
frange = 1./parameters.PeriodRange; 
N_wl = parameters.Wavelengths; % for min. number of wavelengths allowed
Npers = parameters.npers; % Number of periods
xlims = [frange(2) frange(1)]; % limits for plotting
t_vec_all = 1./flip(linspace(frange(1) , frange(2) ,Npers)); % periods at which to extract phase velocity

damp = parameters.damp;
is_normbessel = parameters.is_normbessel;
iswin = parameters.iswin;

minstadist = parameters.minStaDist;
IsFigure = 1;
isfigure2 = 0;
isfigure_snr = 1;

test_funccount_running = 0;
%% Make the initial phase velocity dispersion model

% % calc_Rayleigh_disp
% vec_h = [3 2 4 12]; % Layer thickness
% vec_vs = [1.1 1.2 2.8 3.7 4.6];
% vec_vp = vec_vs.*1.8; vec_pv(1) = 1.5;
% vec_rho = [1.03 1.5 3.02 3.027 3.342];
% vr = mat_disperse(vec_h,vec_rho,vec_vp,vec_vs,1./t_vec_all);
% c = vr(:,1)';
% c_start = c;

% Manually
% c = [3.4837    3.6341    3.7458    3.8223    3.8878    3.9451   4.0013    4.0522    4.0951    4.1337    4.1683    4.2098]; % 'test_1.6win' avg
% 10.0000   11.2063   12.5580   14.0729   15.7704   17.6727   19.8045   22.1935   24.8706   27.8706   31.2325   35.0000
% c = [3    3.2    3.3    3.4464    3.8400    3.9589    4.0097    4.0363    4.0515    4.0600    4.0644    4.0661];

% From MINEOS .q file (https://github.com/jbrussell/MINEOS_synthetics)
qfile = ['./qfiles/Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_ORCAiso_INV.s0to200.q'];
mode = parameters.mode; %wbh put this in params file
if exist('c','var') == 0 % check if phase velocities exist, if not read them in
    [~,~,c_all] = readMINEOS_qfile2(qfile,t_vec_all,mode);
end
%%% wbh try different starting model
c_all = [4.0    4.0    4.0    4.0    4.0    4.0    4.0    4.0]; % try 4 across the board
c_all = [4.5442    4.4953    4.4609    4.4360    4.4174    4.4033    4.3924    4.3837]; % from q file
c_all = [4.5441    4.4953    4.4609    4.4360    4.4174    4.4033    4.3924    4.3837]; % slightly modify
c_all = [4.2    4.2    4.2    4.2    4.2    4.2    4.2    4.2]; % try 4.2 across the board
c_all = [4.2    4.19    4.18    4.17    4.16    4.15    4.14    4.13]; % try decreasing
%c_all = [4.1    4.09    4.08    4.07    4.06    4.05    4.04    4.03]; % try decreasing
c_all = [3.1    3.09    3.08    3.07    3.06    3.05    3.04    3.03]; % try decreasing
c_all = [3.2    3.8    4.0    4.05    3.95    3.9    3.85    3.92]; % approximately J27-J20
%%%
c_start = c_all;
c_all_std = zeros(size(c_all));
%==========================================================%

% LIMITS
if comp{1}(1) == 'R'
    ylims = [3.2 4.5];
elseif comp{1}(1) == 'Z' || comp{1}(1) == 'P'
    ylims = [2.5 5.5];
elseif comp{1}(1) == 'T'
    ylims = [3.5 4.8];
end

% LOAD DATA TO SEE HOW MANY POINTS
%%% --- Load in the ccf --- %%%
        %ccf_path = ['./ccf/',windir,'/fullStack/ccf',comp{1},'/'];
        ccf_path = [parameters.ccfpath,windir,'/fullStack/ccf',comp{1},'/'];
        stalist = parameters.stalist;
        sta1=char(stalist(1,:));
        sta2=char(stalist(2,:));
        sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        if ~exist(filename,'file')
            disp(['not exist ',filename])
        end
        data1 = load(filename);
        npts = length(data1.coh_sum_win);


% input path
%ccf_path = ['./ccf/',windir,'/fullStack/ccf',comp{1},'/'];
ccf_path = [parameters.ccfpath,windir,'/fullStack/ccf',comp{1},'/'];

% output path
%XSP_path = [parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(N_wl),'wl_phv_dir/'];
%%% simple starting model
XSP_path = [parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/',num2str(1/frange(2)),'_',num2str(1/frange(1)),'s_',num2str(N_wl),'wl_phv_dir_Iterate20/'];
%%%

if ~exist(XSP_path)
    if ~exist(parameters.xsppath)
        mkdir(parameters.xsppath);
    end
    if ~exist([parameters.xsppath,windir,'/'])
        mkdir([parameters.xsppath,windir,'/']);
    end
    if ~exist([parameters.xsppath,windir,'/fullStack/'])
        mkdir([parameters.xsppath,windir,'/fullStack/']);
    end
    if ~exist([parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/'])
        mkdir([parameters.xsppath,windir,'/fullStack/Xsp',comp{1},'/']);
    end
    mkdir(XSP_path)
end

% figure output path
if iswin
    %XSP_fig_path = [figDir,windir,'/fullStack/',num2str(N_wl),'wl_phv_dir/TEI19/'];
    XSP_fig_path = [figDir,windir,'/fullStack/',num2str(N_wl),'wl_phv_dir/TEI19_Iterate20/'];
else
    XSP_fig_path = [figDir,windir,'/fullStack/',num2str(N_wl),'wl_phv_dir/TEI19_nowin/'];
end

if ~exist(XSP_fig_path)
    mkdir(XSP_fig_path);
end

warning off; %#ok<WNOFF>

stalist = parameters.stalist;
nsta=parameters.nsta; % number of target stations to calculate for

%%% --- Loop through station 1 --- %%%
for ista1=1:nsta
    
    sta1=char(stalist(ista1,:));
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    
    %%% --- Loop through station 2 --- %%%
    for ista2 = 1: nsta % length(v_sta)
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        % Check to see if we have already done this
        if is_resume && exist([XSP_path,sta1,'_',sta2,'_xsp.mat'])
            disp('Already fit this one!')
            continue
        end
        clear data1 xcorf1 xsp1 filename
        
        %%% --- Load in the ccf --- %%%
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        
        
        data1 = load(filename);
        r1 = distance(data1.stapairsinfo.lats(1),data1.stapairsinfo.lons(1),data1.stapairsinfo.lats(2),data1.stapairsinfo.lons(2),referenceEllipsoid('GRS80'))/1000;
        groupv_max = data1.max_grv;
        groupv_min = data1.min_grv;
        
        if r1 < minstadist
            continue;
        end
        
        % Index wavelength criterion
        I_wl = r1 ./ (t_vec_all .* c_all) > N_wl;
        if sum(I_wl) <= 1
            I_wl(1) = 1;
            I_wl(2) = 1;
        end
        c = c_all(I_wl);
        t_vec = t_vec_all(I_wl);
        c_std = c_all_std(I_wl);
        
        tN = length(t_vec);
        wholesec = npts;
        wvec1 = (2*pi)./t_vec;
        wvec1 = wvec1';
        
        % Get your axis correct
        twloc=1./t_vec;
        twloc = twloc*2*pi;
%         waxis = (frange(1):1/wholesec:frange(2))*2*pi;
        waxis = (1/max(t_vec):1/wholesec:1/min(t_vec))*2*pi;
        
        
        %%% - Get the normalized ccf - %%%
        if iswin
            xcorf1 = data1.coh_sum_win./data1.coh_num;
        else
            xcorf1 = data1.coh_sum./data1.coh_num;
        end
        dumnan = find(isnan(xcorf1)==1);
        
        if length(dumnan) > 10
            disp([sta1,' and ',sta2,'is NaN! Moving on']);
            continue
        end
        
        %N = 10000;
        N = length(xcorf1);
        if length(xcorf1) < N
            disp('Dataset is too short! Moving on')
            continue
        end
        
        xcorf = xcorf1;
        xcorf1 = real(xcorf1(1:N));
        xcorf1(1) = 0;
        
        if isfigure2 
            figure(1)
            T = length(xcorf1);
            dt = 1;
            temp_faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(temp_faxis>0);
            subplot(2,1,1)
            %plot(temp_faxis(ind),smooth(real(xcorf1(ind)),100));
            plot(flip(temp_faxis(ind)),smooth(real(xcorf1(ind)),50));
            xlim([frange(2) frange(1)])
            hold on
            subplot(2,1,2)
            %plot(temp_faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),100),'-r')
            plot(temp_faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),50),'-r');
            xlim([frange(1) frange(2)])
            
        end

        %%% - Convert xcorf into spherical frequency - %%%
        faxis = [0:N-1]*1/wholesec;
        xsp1 = interp1(faxis*2*pi,xcorf1,waxis);

        %xsp1 = smooth(xsp1,50);
        xsp1 = smooth(xsp1,npts_smooth);

        tw1 = ones(1,tN)*r1./c;
        weight  = 1./waxis;
%% WBH Old version        
%         %%% - Invert for the bessel function 2x - %%%
%         options = optimoptions(@lsqnonlin,'TolFun',1e-12,'MaxIter',1500,'MaxFunEvals',1500);
%         
%         %tw2 = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw1],[tw1]*0.8,[tw1]*1.2,options);
%         [tw2,resnorm2,resid2,exitflag2,output2,lambda2,j2] = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw1],[tw1]*0.8,[tw1]*1.2,options);
% %         tw2 = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[],[],options);
%         
%         weight(:) = 1;
%         [tw,resnorm3,res,exitflag3,output3,lambda3,J] = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw2],[tw2]*0.8,[tw2]*1.2,options);
%         disp([tN,length(tw2),length(tw)])
% %         tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[tw2]*0.8,[tw2]*1.2,options);
% %         tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[],[],options);

%% WBH new version to iterate over starting models.
%         options = optimoptions(@lsqnonlin,'TolFun',1e-12,'MaxIter',1500,'MaxFunEvals',1500);
%         c_i = c; %velocity model of iterations... 
%         [tw_final,resnorm_final,resid_final,exitflag_final,output_final,lambda_final,j_final] = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw1],[tw1]*0.8,[tw1]*1.2,options);
%         test_funcCount = 0;
%         % check to see if it stopped bc it hit max fun evals
%         %if output_final.funcCount > options.MaxFunctionEvaluations
%             % we're going to do this until c_i is above some velocity
%             while c_i(1) < 5.0
%                 c_i = c_i + 0.1; % add 0.1 km/s to each phase velocity
%                 tw_new = ones(1,tN)*r1./c_i; %calculate new tw
%                 %do inversion again with new starting model
%                 [tw_i,resnorm_i,resid_i,exitflag_i,output_i,lambda_i,j_i] = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw_new],[tw_new]*0.8,[tw_new]*1.2,options);
%                 % keep these values if they're better than the above
%                 if resnorm_i < resnorm_final
%                     tw_final = tw_i;
%                     resnorm_final = resnorm_i;
%                     resid_final = resid_i;
%                     exitflag_final = exitflag_i;
%                     output_final = output_i;
%                     lambda_final = lambda_i;
%                     j_final = j_i;
%                 end
%                 % let's test to see if the fun count works for breaking the
%                 % loop
%                 if output_i.funcCount < options.MaxFunctionEvaluations
%                     if test_funcCount == 0
%                         %this is the first time
%                         resid_test = resid_final;
%                         test_funcCount = 1;
%                     elseif test_funcCount == 1
%                         %this ought to be every other time...
%                         if resid_test > resid_final
%                             %this would be a problem; a lower resid after
%                             %we would have broken out of the loop
%                             %pause
%                             %don't pause but...
%                             test_funccount_running = test_funccount_running + 1;
%                         end
%                     end
%                     
%                         
%                 end
%             end
%         %end
%% WBH Instead choose a good starting model and add random noise for some set number of iterations

        options = optimoptions(@lsqnonlin,'TolFun',1e-12,'MaxIter',1500,'MaxFunEvals',1500);
        c_i = c; %velocity model of iterations... 
        [tw_final,resnorm_final,resid_final,exitflag_final,output_final,lambda_final,j_final] = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw1],[tw1]*0.8,[tw1]*1.2,options);
        Niter = 20;
        ii = 0;
        tw_starting = tw1;
        dvMagnitude = 0.2; %magnitude of the random perturbation in km/s
        while ii < Niter
            ii = ii + 1;
            % generate random numbers
            dvRand = rand(size(c));
            %center on zero, note now between -0.5 and 0.5
            dvRand = dvRand - 0.5;
            % apply the right magnitude range
            dvRand = dvRand * 2 * dvMagnitude;
            % and add to c
            c_ii = c_i + dvRand;
            % calculate new tw
            tw_new = ones(1,tN)*r1./c_ii; 
            %do inversion again with new starting model
            [tw_i,resnorm_i,resid_i,exitflag_i,output_i,lambda_i,j_i] = lsqnonlin(@(x) besselerr(x,[xsp1],damp,is_normbessel),[tw_new],[tw_new]*0.8,[tw_new]*1.2,options);
            % check to see if this is a better fit...
            if resnorm_i < resnorm_final
                % if it is, store all the new info and re-do
                tw_final = tw_i;
                tw_starting = tw_new;
                resnorm_final = resnorm_i;
                resid_final = resid_i;
                exitflag_final = exitflag_i;
                output_final = output_i;
                lambda_final = lambda_i;
                j_final = j_i;
                % also store the new velocity model
                c_i = c_ii;
            end
        end
        
%% end of bessel fit
        %apply final values
        tw = tw_final;
        res = resid_final;
        J = j_final;
        
        % ESTIMATE ERROR BARS ON MODEL! :: JBR - 2/2020
        % Calculate data variance from residual following Menke QMDA book eq. (4.31)
        % s_d^2 = E / (N-M)
        sigma_d2 = res'*res / (length(res)-length(tw));
        Cov_m = inv(J'*J)*sigma_d2;
        sigma_m_tw = diag(Cov_m).^(1/2);
        % Propogate model error to phase velocity
        % c = r1./tw;  therefore   dc = |r*t^(-2) * dt|
        sigma_m_c = abs(r1.*tw'.^(-2).*sigma_m_tw);
        
        %%% - Set up the variable structure - %%%
        xspinfo.sta1 = sta1;
        xspinfo.sta2 = sta2;
        xspinfo.lat1 = data1.stapairsinfo.lats(1);
        xspinfo.lon1 = data1.stapairsinfo.lons(1);
        xspinfo.lat2 = data1.stapairsinfo.lats(2);
        xspinfo.lon2 = data1.stapairsinfo.lons(2);
        
        xspinfo.r = r1;
        xspinfo.tw = tw;
        xspinfo.xsp = xsp1;
        xspinfo.xsp_norm = xsp1./abs(hilbert(xsp1));
        xspinfo.coherenum = data1.coh_num;
        err = besselerr(tw,xsp1,damp,is_normbessel);
        err = err(1:length(waxis));
        if is_normbessel
            xspinfo.sumerr = sum(err.^2)./sum((xspinfo.xsp_norm./weight(:)).^2);
        else
            xspinfo.sumerr = sum(err.^2)./sum((xsp1./weight(:)).^2);
        end
        xspinfo.err = err./weight(:);
        xspinfo.tw1 = tw1;
        xspinfo.twloc = twloc;
        xspinfo.c = r1./tw;
        xspinfo.c_std = sigma_m_c;
        xspinfo.per = 1./(twloc/2/pi);
        xspinfo.c_start = c_start;
        xspinfo.c_std_start = c_all_std;
        xspinfo.per_start = t_vec_all;
        xspinfo.isgood_wl = I_wl;
        
        data = r1./tw;
        

        %% %%% Calculate SNR %%%
        xcorf1 = data1.coh_sum./data1.coh_num;
        xcorf1_filtered = tukey_filt( xcorf1,[min(t_vec) max(t_vec)],1,0.25 );
        [snr, signal_ind] = calc_SNR(xcorf1_filtered,groupv_min,groupv_max,r1,isfigure_snr);
        %%

        xspinfo.filename = filename;
        xspinfo.snr = snr;
        
        % Calculate the predicted bessel function from the initial model
        disp([filename,' fitted'])
        if IsFigure
            if 0 % plot initial bessel
                tw_init = interp1(twloc,tw1(1:tN),waxis,'linear');
                x_init = waxis.*tw_init;
                A = 1;
                binit = besselj(0,x_init)*A;
                binit = binit./mean(abs(binit)).*mean([abs(xsp1)]);
                plot(waxis/2/pi,binit,'-k','linewidth',2); hold on;
            end
            
            f3 = figure(3); clf; hold on; 
            set(gcf,'color','w','Position',[289     1   517   704]);

            ax1 = subplot(3,1,1);
            plot_SNR(xcorf1_filtered,groupv_min,groupv_max,r1,ax1);
            set(ax1,'box','off');

            % REAL PART (J0)
            subplot(3,1,2); box off;
            tww = interp1(twloc,tw(1:tN),waxis,'linear');
            x = waxis.*tww;
            A = 1;
            b = besselj(0,x)*A;
            b = b./mean(abs(b)).*mean([abs(xsp1)]);           
            T = length(data1.coh_sum);
            dt = 1;
            faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(faxis>0);
            plot(faxis(ind),smooth(real(data1.coh_sum_win(ind)/data1.coh_num),npts_smooth),'-k','linewidth',3); hold on;
            if ~iswin
                plot(waxis/2/pi,xsp1,'-b','linewidth',1);
            end
            plot(waxis/2/pi,b,'-r','linewidth',2); hold on; 
            xlim(xlims);
            xlims1 = get(gca,'XLim');
            ylabel('J_{0}','fontsize',16);
            ax1 = get(gca);
            dx = 1;
            dy = 0.95;
            text(gca,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
                sprintf('%.4f',xspinfo.sumerr),'color',[0 0 0],'fontsize',14);
            set(gca,'fontsize',16,'linewidth',1.5);
            box off;
            title(['Distance : ',num2str(r1),'km'],'fontsize',16);
                
            hold on
            subplot(3,1,3);
            %errorbar(twloc/2/pi,r1./tw1,c_std,'o-','color',[0.5 0.5 0.5],'linewidth',2);hold on;
            errorbar(twloc/2/pi,r1./tw_starting,c_std,'o-','color',[0.5 0.5 0.5],'linewidth',2);hold on;
            errorbar(twloc/2/pi,r1./tw,sigma_m_c*2,'o-','color',[1 0 0 ],'linewidth',2);
%             errorbar(twloc/2/pi,r1./tw,xspinfo.err,'ro-','linewidth',2);
            title([sta1,'-',sta2],'fontsize',16)
            xlabel('Frequency (Hz)','fontsize',16);
            ylabel('Phase Velocity (km/s)','fontsize',16);
            set(gca,'fontsize',16,'linewidth',1.5);
            ylim(ylims);
            xlim(xlims1);
            box off;
            
            % Plot normalized bessel functions
            if is_normbessel
                figure(4); clf;
                b_dat = smooth(real(data1.coh_sum_win(ind)/data1.coh_num),npts_smooth);
                plot(faxis(ind),b_dat./abs(hilbert(b_dat)),'k','linewidth',3); hold on;
                plot(waxis/2/pi,b./SmoothAnalyticEnv(waxis/2/pi,b),'-r','linewidth',2); hold on;
                xlim(xlims);
                xlims1 = get(gca,'XLim');
                ylabel('J_{0}','fontsize',16);
                ylim([-2 2]);
            end
            
            if isfigure2
            f12 = figure(12);
            clf
            T = length(data1.coh_sum);
            dt = 1;
            faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(faxis>0);
            plot(faxis(ind),smooth(real(data1.coh_sum(ind)/data1.coh_num),npts_smooth),'-r');
            xlim([frange(1) frange(2)])
            end
            psfile = [XSP_fig_path,'Xsp_',comp{1}(1),'_',sta1,'_',sta2,'_J0J1.pdf'];
            %print('-dpsc2',psfile);
            drawnow
            if isoutput
                %save2pdf(psfile,f3,250);
                saveas(f3,psfile);
            end

            
%             pause;
        end
        if isoutput
            save(sprintf('%s/%s_%s_xsp.mat',XSP_path,sta1,sta2),'xspinfo','twloc','waxis');
        end
        
%         pause;
        
        
    end %end of station j
end  %end of station i