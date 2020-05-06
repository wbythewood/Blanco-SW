function spitplots_correctevent(dt,T1,T2,Zraw,H1raw,H2raw,Praw,taxis,eventid,netsta,corrseis,f,NFFT)

fn = 1/2/dt;
[b,a]=butter(2,[1/fn/T2,1/fn/T1]);
Z_filt  = filtfilt(b,a,Zraw);
H1_filt  = filtfilt(b,a,H1raw);
H2_filt  = filtfilt(b,a,H2raw);
P_filt  = filtfilt(b,a,Praw);

figure(101);clf
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
subplot(411)
%plot(taxis,Z_filt,'-k');
plot(taxis,Zraw,'-k');
xlim([min(taxis),max(taxis)]);
title(sprintf('%s %s Z',netsta, eventid));
subplot(412)
%plot(taxis,H1_filt,'-k');
plot(taxis,H1raw,'-k');
xlim([min(taxis),max(taxis)]);
title(sprintf('%s H1',eventid));
subplot(413)
%plot(taxis,H2_filt,'-k');
plot(taxis,H2raw,'-k');
xlim([min(taxis),max(taxis)]);
title(sprintf('%s H2',eventid));
subplot(414)
%plot(taxis,P_filt,'-k');
plot(taxis,Praw,'-k');
title(sprintf('%s P',eventid));
xlim([min(taxis),max(taxis)]);

figure(102); clf
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
subplot(length(corrseis),1,length(corrseis))
plot(taxis,Z_filt,'-k');
xlim([min(taxis),max(taxis)]);
title(sprintf('Original Z'));
for is = 1:length(corrseis)
    if corrseis(is).isgood ==1
        col = [0 0 0];
    else
        col = [1 0 0];
    end
    if strcmp(corrseis(is).label(1),'Z')==1
        seiscmp = Z_filt;
    elseif strcmp(corrseis(is).label(1),'1')==1
        seiscmp = H1_filt;
    elseif strcmp(corrseis(is).label(1),'2')==1
        seiscmp = H2_filt;
    elseif strcmp(corrseis(is).label(1),'P')==1
        seiscmp = P_filt;
    end
    subplot(length(corrseis),1,is)
%     plot(taxis,seiscmp,'-','LineWidth',.5,'Color',[.5 .5 .5]); hold on
    plot(taxis,filtfilt(b,a,corrseis(is).timeseries),'-','LineWidth',.5,'Color',col)
    xlim([min(taxis),max(taxis)]);
    title(sprintf('%s %s %s',netsta,eventid,corrseis(is).label));
end

% wbh add figure to see spectral properties of corrected traces
figure(103); clf
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
subplot(length(corrseis),1,length(corrseis))
for is = 1:length(corrseis)
    subplot(length(corrseis),1,is)
    psd = disp2accel_psd(corrseis(is).timeseries,NFFT,dt);
    semilogx(f,10*log10(psd),'-k','linewidth',0.5); 
    title(sprintf('%s %s %s',netsta,eventid,corrseis(is).label))
    grid on; set(gca,'layer','top');
    xlim([1/250,fn]);
    xlabel('Freq (Hz)');
    ylim([-200,-80]);
    ylabel('Power (DB)');
end

% wbh add figure to see spectra and traces next to each other
figure(104); clf
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
subplot(length(corrseis),2,length(corrseis))
for is = 1:length(corrseis)
    if corrseis(is).isgood == 1
        col = [0 0 0];
    else
        col = [1 0 0];
    end
    subplot(length(corrseis),2,2*is)
    psd = disp2accel_psd(corrseis(is).timeseries,NFFT,dt);
    semilogx(f,10*log10(psd),'-k','linewidth',0.5); 
    title(sprintf('%s %s %s',netsta,eventid,corrseis(is).label))
    grid on; set(gca,'layer','top');
    xlim([1/250,fn]);
    xlabel('Freq (Hz)');
    ylim([-200,-80]);
    ylabel('Power (DB)');
    subplot(length(corrseis),2,2*is-1)
    plot(taxis,filtfilt(b,a,corrseis(is).timeseries),'-','LineWidth',.5,'Color',col)
    %plot(taxis,corrseis(is).timeseries,'-','LineWidth',.5,'Color',col)
    xlim([min(taxis),max(taxis)]);
    xlabel('Time (s)')
    title(sprintf('%s %s %s',netsta,eventid,corrseis(is).label));
end

return