function [max_coh_fine,max_or_fine] = spec_orient(spectrum_Z,spectrum_H1,spectrum_H2,cspectrum_Z,hangs,tiltfreq,f,isfig,dayid,isgoodwin,NFFT,dt)

c = colormap('jet');

hangint = hangs(2)-hangs(1);

days = dayofyear(str2num(dayid(1:4)),str2num(dayid(5:6)),str2num(dayid(7:8)),0,0);
cc=interp1(1:64,c,((days)/(365))*63+1);

ph_points = zeros(1,2);
for ih = 1:length(hangs)
    hang = hangs(ih);
    cang = cos(hang*pi/180);
    sang = sin(hang*pi/180);

    chh_stack=zeros(1,length(f))';
    czz_stack=zeros(1,length(f))';
    chz_stack=zeros(1,length(f))';
    nwin_stack = 0;
    for ista = 1:length(isgoodwin)
        if isgoodwin(ista)==0
            continue
        end
        spec_H = sang.*spectrum_H2(ista,:)'+cang.*spectrum_H1(ista,:)'; %rotated horizontal spectra and cross-spectra

        chh = abs(spec_H).^2*2/(NFFT*dt);
        czz = abs(spectrum_Z(ista,:)').^2*2/(NFFT*dt);
        chz = spec_H.*(cspectrum_Z(ista,:)')*2/(NFFT*dt);
        % add a matrix that keeps track of all the cross spectra for each
        % time window
        % (freq elements in spectrum, angle from H1, window number)
        chz_matrix(:,ih,ista) = chz;

        % stack the cross spectra
        chh_stack=chh_stack+chh;
        czz_stack=czz_stack+czz;
        chz_stack=chz_stack+chz;
        nwin_stack = nwin_stack +1;
    end

    %Normalization of the stack
    chh_stack = chh_stack/nwin_stack;
    czz_stack = czz_stack/nwin_stack;
    chz_stack = chz_stack/nwin_stack;

    % Coherence
    coh_stack = abs(chz_stack).^2./(chh_stack.*czz_stack);
    % Phase *of the daily stack, i.e., average*
    ph_stack = 180/pi.*atan2(imag(chz_stack),real(chz_stack));
    % make the phase [0,360] instead of [-180,180]
    %ph_stack(ph_stack < 0) = ph_stack(ph_stack<0) + 360;
    % Admittance
    ad_stack = abs(chz_stack)./chh_stack;

    % plots for debugging, make sure rotation is operating correctly
%     figure(95)
%     subplot(211)
%     semilogx(f,(smooth(coh_stack,40)+ih),'-k');
%     hold on
%     semilogx(f,ih*ones(size(f)),'-r');
%     xlim([0.005 0.035])
%
%     figure(95)
%     subplot(212)
%     semilogx(f,(ph_stack+ih*180),'.k');
%     hold on
%     semilogx(f,180*ih*ones(size(f)),'-r');
%     xlim([0.005 0.035])
    % end debugging plots

    % looking for average coherence value within frequency range

    [fmin,idx_flo] = min(abs(f-tiltfreq(1)));
    [fmin,idx_fhi] = min(abs(f-tiltfreq(2)));
    av_coh_coarse(ih) = mean(coh_stack(idx_flo:idx_fhi));
    av_ph_coarse_old(ih) = mean(abs(ph_stack(idx_flo:idx_fhi)));
    % now that it's all positive this shouldn't matter... but before when
    % phase went -180 to 180, mean of abs gave answers centered around 90
    %av_ph_coarse(ih) = abs(mean(ph_stack(idx_flo:idx_fhi)));
    av_ph_coarse(ih) = mean(ph_stack(idx_flo:idx_fhi));

end

% for every corss spectrum, find phase...
chz_matrix = 180/pi.*atan2(imag(chz_matrix),real(chz_matrix));
% only look within specific frequency band...
chz_matrix = chz_matrix(idx_flo:idx_fhi,:,:);

% choose a single random window...
ph_matrix = chz_matrix(:,:,3);
% or average across time windows... may not be sensible
%ph_matrix = mean(chz_matrix,3);
% choose a random single frequecy...
ph_matrix = ph_matrix(8,:);
% or average over all frequencies... but this may not give sensible answers
%ph_matrix = mean(ph_matrix,1);

% average through one dimension... 1=freq; 2=hang, 3=window
reduced_ph_matrix = mean(chz_matrix,3);
reduced_ph_matrix2 = squeeze(mean(chz_matrix,1));
% ...and choose a random value in the other
reduced_ph_matrix = reduced_ph_matrix(8,:);
reduced_ph_matrix2 = reduced_ph_matrix2(:,3).';

% save dimensions of matrix...
NumWindows = size(chz_matrix,3);
NumFreq = size(chz_matrix,1);

PhasePoints = zeros(1,2);
% vector with [hang,phase] for heat map
% for iWin = 1:NumWindows
%     for iFreq = 1:NumFreq
%         for ih = 1:length(hangs)
%             iPhase = [chz_matrix(iFreq,ih,iWin),hangs(ih)];
%             PhasePoints = [PhasePoints; iPhase];
%         end
%     end
% end


% plotting the dependence of phase and coherence on angle
if isfig ==1
    % color scaled by day of year
    figure(105)

    subplot(1,2,1); hold on
    plot(hangs,av_coh_coarse,'LineWidth',.5,'Color',cc);
    title(sprintf('Coherence, %.3f - %.3f Hz', tiltfreq(1),tiltfreq(2)));
    xlabel('Angle from H1'); ylabel ('Coherence'); ylim([0 1]); xlim([0 360])

    subplot(1,2,2); hold on
    plot(hangs,av_ph_coarse_old,'LineWidth',.5,'Color',cc);
    title(sprintf('Average Phase, %.3f - %.3f Hz', tiltfreq(1),tiltfreq(2)));
    xlabel('Angle from H1'); ylabel ('Phase'); ylim([-200 200]); xlim([0 360])
end

% wbh understanding what's going on
figure(106)

subplot(2,2,1); hold on
plot(hangs,av_ph_coarse_old,'LineWidth',.5,'Color',cc);
title(sprintf('Old Average Phase, mean(abs())'));
xlabel('Angle from H1'); ylabel ('Phase'); ylim([-200 200]); xlim([0 360]);

subplot(2,2,2); hold on
%plot(hangs,reduced_ph_matrix,'LineWidth',.5,'Color',cc);
%title(sprintf('Average through t windows, choose f'));
plot(hangs,reduced_ph_matrix2,'LineWidth',.5,'Color',cc);
title(sprintf('Average through f, choose t window'));
xlabel('Angle from H1'); ylabel ('Phase'); ylim([-200 200]); xlim([0,360])

%wbh

subplot(2,2,3); hold on
title(sprintf('Individual phase measurements'));
xlabel('Angle from H1'); ylabel ('Phase'); ylim([-200 200]); xlim([0,360]);

subplot(2,2,4); hold on
title(sprintf('abs(individual phase measurements)'));
xlabel('Angle from H1'); ylabel ('Phase'); ylim([-200 200]); xlim([0,360]);

% lol this prints way too many plots
for iWin = 1:NumWindows
    subplot(2,2,3);
    plot(hangs,chz_matrix(5,:,iWin),'LineWidth',.5,'Color',cc);
    subplot(2,2,4);
    plot(hangs,abs(chz_matrix(5,:,iWin)),'LineWidth',.5,'Color',cc);
end
% try something else... 
    
%plot(hangs,ph_matrix,'LineWidth',.5,'Color',cc);

 


% figure(107)
% hold on;
% scatter(PhasePoints(:,1),PhasePoints(:,2))

idxph = find(abs(av_ph_coarse_old)>90);

[max_coh_coarse,idx]=max(av_coh_coarse(idxph),[],2);
max_or_coarse = hangs(idxph(idx));

% fine grid search to get best orientation
hangs2 = [max_or_coarse-hangint+1:max_or_coarse+hangint-1];
for ih = 1:length(hangs2)
    hang = hangs2(ih);
    if hang == max_or_coarse;
        continue
    end
    cang = cos(hang*pi/180);
    sang = sin(hang*pi/180);

    % Initialize output structures
    chh_stack=zeros(1,length(f))';
    czz_stack=zeros(1,length(f))';
    chz_stack=zeros(1,length(f))';
    nwin_stack = 0;
    for ista = 1:length(isgoodwin)
        if isgoodwin==0
            continue
        end
        spec_H = sang.*spectrum_H2(ista,:)'+cang.*spectrum_H1(ista,:)'; %rotated horizontal spectra and cross-spectra

        chh = abs(spec_H).^2*2/(NFFT*dt);
        czz = abs(spectrum_Z(ista,:)').^2*2/(NFFT*dt);
        chz = spec_H.*(cspectrum_Z(ista,:)')*2/(NFFT*dt);

        chh_stack=chh_stack+chh;
        czz_stack=czz_stack+czz;
        chz_stack=chz_stack+chz;
        nwin_stack = nwin_stack +1;
    end

    %Normalization
    chh_stack = chh_stack/nwin_stack;
    czz_stack = czz_stack/nwin_stack;
    chz_stack = chz_stack/nwin_stack;

    % Coherence
    coh_stack = abs(chz_stack).^2./(chh_stack.*czz_stack);
    % Phase
    ph_stack = 180/pi.*atan2(imag(chz_stack),real(chz_stack));
    % Admittance
    ad_stack = abs(chz_stack)./chh_stack;


    [fmin,idx_flo] = min(abs(f-tiltfreq(1)));
    [fmin,idx_fhi] = min(abs(f-tiltfreq(2)));
    av_coh(ih) = mean(coh_stack(idx_flo:idx_fhi));
    %av_ph(ih) = mean(abs(ph_stack(idx_flo:idx_fhi)));
    av_ph(ih) = mean(ph_stack(idx_flo:idx_fhi));
end

[max_coh_fine,idx]=max(av_coh,[],2);
max_or_fine = hangs2(idx);

max_coh_fine = max_coh_fine;
max_or_fine = max_or_fine;

if isfig ==1
figure(90)
    subplot(211); hold on
    plot(days,max_or_fine,'o','MarkerFaceColor',cc,'MarkerSize',5,'MarkerEdgeColor','none');
    title('Orientation of Maximum Coherence'); xlabel('Days Since January 1'); ylabel('Degrees from H1'); ylim([0 360])
    subplot(212); hold on
    plot(days,max_coh_fine,'o','MarkerFaceColor',cc,'MarkerSize',5,'MarkerEdgeColor','none');
    title('Value of Maximum Coherence'); xlabel('Days Since January 1'); ylabel('Coherence'); ylim([0 1])
end

return
