% JBR - 1/9/16
%
% Plot phase delay (dtp) as a function of epicentral distance difference
% (ddist) similar to Figure 6 in the ASWMS manual.
%

clear;
addpath('../')

setup_parameters;

periods = parameters.periods;

is_fig = 1;
% --------------------------

eventfile = parameters.eventfile;
workingdir = parameters.workingdir;
aswmsdir = parameters.ASWMSDir;
evmat = [workingdir,'data/eventmat/'];
CSmat = [aswmsdir,'CSmeasure/'];
fig_PATH = [workingdir,'figures/'];

maxstadist = parameters.maxstadist;
xlims = [-maxstadist maxstadist];

if ~exist(fig_PATH)
    mkdir(fig_PATH);
end
%% LOAD DATA STRUCTURES

% LOAD MEASUREMENTS .mat
obs_CSpath = [CSmat];
obs_CSfiles = dir([obs_CSpath,'*.mat']);
numCS = length(obs_CSfiles); % one for each event

for iev = 1:numCS
    load([obs_CSpath,obs_CSfiles(iev).name]); %loads structure called "eventcs"
    events_obs(iev).eventcs = eventcs;
end

%% Plot phase delays

for ev = 1:numCS
%Loop over periods for a given event
max_dtp = -1000;
min_dtp = 1000;

% Data
fig58 = figure(58); clf; hold on; set(gcf, 'Color', 'w'); box on;
numper = length(periods);
clr_per = jet(numper);
lgd = {};
h = zeros(numper,1);
for iper = 1:numper
    useper = numper+1-iper;
    useper = iper;
    num_measures = length(events_obs(ev).eventcs.CS);
    for imeas = 1:num_measures
        if events_obs(ev).eventcs.CS(imeas).isgood(useper) ~= 1 %QC
        %if events_obs(ev).eventcs.CS(imeas).isgood(iper) == -2 %QC high tp error
            plot(events_obs(ev).eventcs.CS(imeas).ddist,events_obs(ev).eventcs.CS(imeas).fiterr(useper),'o','color',[.8 .8 .8],'linewidth',1);
        else
            if events_obs(ev).eventcs.CS(imeas).fiterr(useper) < min_dtp
                min_dtp = events_obs(ev).eventcs.CS(imeas).fiterr(useper);
            end
            if events_obs(ev).eventcs.CS(imeas).fiterr(useper) > max_dtp
                max_dtp = events_obs(ev).eventcs.CS(imeas).fiterr(useper);
            end
            h(useper) = plot(events_obs(ev).eventcs.CS(imeas).ddist,events_obs(ev).eventcs.CS(imeas).fiterr(iper),'+','color',clr_per(iper,:),'linewidth',1);
        end
    end
    lgd{useper} = [num2str(periods(useper)),'s'];
end
xlabel('Epicentral Distance Difference (km)','fontsize',12);
xlim(xlims);
ylabel('Fit error','fontsize',12);
title([num2str(events_obs(ev).eventcs.id),' Phase Delays'],'fontsize',12);
legend(h,lgd,'location','northeastoutside');
ax = gca;
set(ax,'ylim',[min_dtp max_dtp],'linewidth',2,'fontsize',16);
drawnow;
% pause;

%% EXPORT FIGURES
if is_fig == 1
%     export_fig(fig55,[fig_PATH,'obs_',num2str(events_obs(ev).eventcs.id)],'-pdf','-painters');
    %save2pdf([fig_PATH,'obs_',num2str(events_obs(ev).eventcs.id)],fig55)
    ofn = [fig_PATH,'obs_rev_',num2str(events_obs(ev).eventcs.id),'.pdf'];
    saveas(fig58,ofn)
%     export_fig(fig66,[fig_PATH,'synth_',num2str(events_synth(ev).event.id)],'-pdf','-painters');
end
%pause;
end