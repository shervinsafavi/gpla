%% basics

% clear all
% ignit_gpla
pds = ignit();

% load(fullfile(dirPath_gpla.data, 'handies', 'hpcSmltHandyVars_wAnalyticLfp.mat'))
load(fullfile(pds.ldat, 'handies', 'hpcSmltHandyVars_wAnalyticLfp.mat'))
% load(fullfile(dirPath_gpla.data, 'dataHPsim', 'Session(full)_10_CA3CA1Network.mat'))
load(fullfile(pds.ldat, 'dataHPsim', 'Session(full)_10_CA3CA1Network.mat'))


sc_threshold = 500;
nfreq = numel(freqBands);
caseName = 'wrc_rmth';

%%
sameElecCheckInfo = [];
statTestInfo = []; % no stat test will be applied

%% GPLA
clear lfpVec spkVec
iSV = 1; % this will change in later part of the analysis

for iFreq = 1 : nfreq
    lfpPh = lfpAnalSig_allSepFreq{iFreq};

    [lfpVec.(caseName)(:,iFreq), spkVec.(caseName)(:,iFreq), gPLV.(caseName)(iFreq), ...
     ~, ~, rawSvdStuff.(caseName)(iFreq)] = ...
        tngpla(spkTrains, lfpPh ,[], sc_threshold, [], [] , 1, ...
               statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);
end

%% ~ plot
cell2mat(freqBands')

cellCol = [...
    1 0 0;...
    0 0 1; ...
    .8 .5 0; ...
    0 .5 .8 ...
          ];

nType = 4;

% nR = nfreq;
% nC = nType;

cellLabels = {'CA3 pyr','CA3 int','CA1 pyr','CA1 int'}

isv = 1;
iFreq = 3;
% for iFreq = 1 : nfreq
tsv = spkVec.(caseName)(:,iFreq, isv);

for kType = 1 : nType
    polarplot(angle(tsv(Sig.cell_label == kType)), abs(tsv(Sig.cell_label == kType)), ...
              '.', 'color',cellCol(kType,:), 'markersize',12);
    hold on
end

thetaticks(0:45:315)
% rRange = get(gca, 'rlim')
% tn = find_closestEvenNum1(rRange(2));
%
tn = .18;
rticks([0 tn/3 2*tn/3 tn]);
rlim([0 tn])
ax = gca;

rruler = ax.RAxis;
% rruler.Label.String = rlabel;

ax.RAxisLocation = 67.5000;
% ax.RAxisLocation = 22.5;

ax.GridAlpha = .2;
ax.LineWidth = 1.5;


%%
lh = legend('CA3 pyr','CA3 int','CA1 pyr','CA1 int','location','east')
displaceFigureStuff(lh, [.155 NaN NaN NaN])
set(lh, 'box', 'off')
set(gca, 'fontsize', 8.5)

