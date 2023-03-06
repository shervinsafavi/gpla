% todo:
% 1- name the PLV similar values in other figs
% 2- values shpuld be normzilied

clear all
% ignit_gpla

run('../vizConventions')

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
clf
nR = 9;
nC = 4;

for iFreq = 1 : 8
    % iFreq = 3;
    % subplot2(nR,nC, [2 3], [1 2]);
    subplot2(nR,nC, iFreq+1, 2+1);
    image_mapArray(lfpVec.uate2010(:,iFreq),utahMaps(iDataset).LFP.map);
    subplot2(nR,nC, iFreq+1, 3+1);
    image_mapArray(spkVec.uate2010(:,iFreq),utahMaps(iDataset).MU.map);

end

for iFreq = 1 : 8
    subplot2(nR,nC, iFreq+1, 2);
    polarplot(angle(lfpVec.uate2010(:,iFreq)), ...
              abs(lfpVec.uate2010(:,iFreq)), 'k.', 'MarkerSize', 5, 'color',vc.lfpv); 
    hold on
    polarplot(angle(spkVec.uate2010(:,iFreq)), ...
              abs(spkVec.uate2010(:,iFreq)), 'k.', 'MarkerSize', 5, 'color',vc.spkv); 
end

theta = 0 : .01 : 2*pi;

for iFreq = 1 : 8
    svec = spkVec.uate2010(:,iFreq);
    [mu, kappa] = getMacroPic(svec);
    y = vmpdf1(theta, mu, kappa);

    subplot2(nR,nC, iFreq+1, 1);
    polarplot(theta, y * gPLV.uate2010(iFreq), ...
              'color',vc.spkv);

    svec = lfpVec.uate2010(:,iFreq);
    [mu, kappa] = getMacroPic(svec);
    y = vmpdf1(theta, mu, kappa);

    % subplot2(nR,nC, iFreq+1, 1); %
    hold on
    polarplot(theta, y * gPLV.uate2010(iFreq), ...
              'color',vc.lfpv);


    % image_mapArray(lfpVec.uate2010(:,iFreq),utahMaps(iDataset).LFP.map);
    % image_mapArray(spkVec.uate2010(:,iFreq),utahMaps(iDataset).MU.map);
end

% hlp_savefig('~/ownCloud/research/nnr/figures_nnr/toHave_20180618/fig5_UATE2010');

%%
fsp = '.'
svflg = 1;
fn = 'fig5_an1_UATE2010res';
spf(fsp, fn, svflg)

