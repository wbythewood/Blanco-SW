function plot_cleanstaspectra(spect,coh_stack,ph_stack,ad_stack,cc,station,f,maxpow,minpow);

comporder = {'Z','H1','H2','P'};
plotorder = {'1Z','2Z','PZ','12','1P','2P','12'};

% Plotting Power Spectra
figure(1)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);
for ip=1:4
    subplot(4,1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    loglog(f,spect(:,:,ip),'-','LineWidth',.5)%,'Color',cc')
%     hold on
    title(sprintf('%s-component, Station: %s',comporder{ip},station));
    xlabel('Frequency (Hz)')
    ylabel('Power (db)')
    set(gca,'yscale','log','xscale','log');
    xlim([10^-4 max(f)]); 
    if isnan(minpow(ip))
        ylim([ 0 1])
    elseif minpow(ip)==maxpow(ip)
        ylim([ 0 1])
    else
    ylim([minpow(ip) maxpow(ip)]);
    end
    box on
end

figure(2)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

figure(3)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

figure(4)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 10.5]);

for ip=1:6
    figure(2)
    subplot(6,1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    semilogx(f',coh_stack(:,:,ip),'-','LineWidth',.5);%,'Color',cc);
    hold on
    title(sprintf('%s Coherence: %s',station, plotorder{ip}))
    xlabel('Frequency (Hz)')
    ylabel('Coherence')
    ylim([0 1]); xlim([10^-4 max(f)]);
    set(gca,'xscale','log');
    box on
    
    figure(3)
    subplot(6,1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    semilogx(f',ph_stack(:,:,ip),'o','MarkerSize',1);
    hold on
    title(sprintf('%s Phase: %s',station,plotorder{ip}))
    xlabel('Frequency (Hz)')
    ylabel('Phase')
    xlim([10^-4 max(f)]);
    set(gca,'xscale','log');
    box on
    
    figure(4)
    subplot(6,1,ip)
    set(gca,'ColorOrder',cc,'NextPlot','replacechildren');
    loglog(f',ad_stack(:,:,ip),'-','LineWidth',.5)
    hold on
    title(sprintf('%s Admittance: %s',station,plotorder{ip}))
    xlabel('Frequency (Hz)')
    ylabel('Admittance')
    xlim([10^-4 max(f)]);
    set(gca,'yscale','log','xscale','log');
    box on
    
end

return