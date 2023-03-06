% todo:
% 1- name the PLV similar values in other figs
% 2- values shpuld be normzilied

clear all
% ignit_gpla

vc = vizConventions()

%%
freqBands = [...
    3 5;...
    5 15; ...
    15 30; ...
    30 50; ...
    50 70; ...
    70 105; ...
    105 125; ...
    125 140, ...
    ];
freqBandsLab = {...
    '3-5';...
    '5-15'; ...
    '15-30'; ...
    '30-50'; ...
    '50-70'; ...
    '70-105'; ...
    '105-125'; ...
    '125-140', ...
    };

iDataset = 1

%% Figure 5: Utah array data
% This figure depict the application of the method on data recorded from
% monkey PFC

%%
% Prepare the processed data: 
load([dirPath_gpla.data, filesep, 'uate2010_jtrBsdsigTest_20180504.mat'])

% raw data: 
load([dirPath_gpla.data, filesep, 'metaData_UA_20161018.mat'])
load([dirPath_gpla.data, filesep, 'spikeData_UA_20161018.mat'])
load([dirPath_gpla.data, filesep, 'lfpData_UA_20161019.mat'])

% tmp stuff
% lfpVec.waveIllu{iCase}{iBeta}, spkVec.waveIllu{iCase}{iBeta}, gPLV.waveIllu{iCase}{iBeta}, ...
%             ~, gPLV_stats.waveIllu{iCase}{iBeta}, ...

tmp = gPLV; clear gPLV
gPLV.uate2010 = tmp; clear tmp

tmp = lfpVec; clear lfpVec  
lfpVec.uate2010 = tmp; clear tmp

tmp = spkVec; clear spkVec
spkVec.uate2010 = tmp; clear tmp

tmp = gPLV_stats; clear gPLV_stats
gPLV_stats.uate2010 = tmp; clear tmp

tmp = PLV_stats; clear PLV_stats
PLV_stats.uate2010 = tmp; clear tmp

%% 

iFreq = 3;

nR = 3;
nC = 4;

subplot2(nR,nC, 1, [1 2]);
bar(gPLV.uate2010 ,'k');
pbaspect([1 .48 1])
% plot_sigStar(abs(gPLV.(tmpCh)), pVal_gpla.(tmpCh))
% hold all
% plot(6, max(abs(gPLV.(tmpCh)))+.35, 'rv', 'MarkerSize',20, 'MarkerFaceColor','r')
% 
% % ylim([0 4]);
% xlim([0 9])
set(gca,'xticklabel',freqBandsLab)
xtickangle(45)
xlabel('Frequency band (Hz)')
ylabel('gPLV')
% set(gca, 'fontSize', tmpFontSize)
grid on

subplot2(nR,nC, 1, 3);
polarplot(angle(lfpVec.uate2010(:,iFreq)), ...
    abs(lfpVec.uate2010(:,iFreq)), 'k.', 'MarkerSize', 5, 'color',vc.lfpv); 
subplot2(nR,nC, 1, 4);
polarplot(angle(spkVec.uate2010(:,iFreq)), ...
    abs(spkVec.uate2010(:,iFreq)), 'k.', 'MarkerSize', 5, 'color',vc.spkv); 

subplot2(nR,nC, [2 3], [1 2]);
image_mapArray(lfpVec.uate2010(:,iFreq),utahMaps(iDataset).LFP.map)
subplot2(nR,nC, [2 3], [3 4]);
image_mapArray(spkVec.uate2010(:,iFreq),utahMaps(iDataset).MU.map)

% hlp_savefig('~/ownCloud/research/nnr/figures_nnr/toHave_20180618/fig5_UATE2010');
