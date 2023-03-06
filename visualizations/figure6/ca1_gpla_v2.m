% main difference with v1 is the location of subplots

%%
load(fullfile(dirPath_gpla.data, 'handies', 'hpcSmltHandyVars_wAnalyticLfp.mat'))
load(fullfile(dirPath_gpla.data, 'dataHPsim', 'Session(full)_10_CA3CA1Network.mat'))

sc_threshold = 500;
nfreq = numel(freqBands);

%%
sameElecCheckInfo = [];
statTestInfo = []; % no stat test will be applied

%%
cellLabels = {'CA3 pyr','CA3 int','CA1 pyr','CA1 int'}
cellCol = [...
    1 0 0;...
    0 0 1; ...
    .8 .5 0; ...
    0 .5 .8 ...
          ];

%% pick CA1
regionLfpsInd = 1 : 32;
selCellID = [3 4];

regionSpkInd = find(Sig.cell_label == selCellID(1) | Sig.cell_label == selCellID(2));

% LFP 
lfpAnalSig_allSepFreq_old = lfpAnalSig_allSepFreq;
clear lfpAnalSig_allSepFreq
for ifreq = 1 : nfreq
    lfpAnalSig_allSepFreq{ifreq} = lfpAnalSig_allSepFreq_old{ifreq}(regionLfpsInd, :,:);
end
n.nCh = 32;

clear spkTrains_old lfpAnalSig_allSepFreq_old

%% GPLA
caseName = 'ca1anal';
% clear lfpVec spkVec
iSV = 1; % this will change in later part of the analysis

for iFreq = 1 : nfreq
    lfpPh = lfpAnalSig_allSepFreq{iFreq};

    [lfpVec.(caseName)(:,iFreq), spkVec.(caseName)(:,iFreq), gPLV.(caseName)(iFreq), ...
     ~, ~, rawSvdStuff.(caseName)(iFreq)] = ...
        tngpla(spkTrains, lfpPh ,[], sc_threshold, regionSpkInd, [] , 1, ...
               statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);
end

%% plot gPLV
% subplot2d(nR,nC, sdi+1, [1 2]);
k = 3;
subplot2d(nR,nC, k, [3 4]-2);
bar(gPLV.(caseName), 'k')

%% ~ plot spike vectors
cell2mat(freqBands')

nType = 4;

cellLabels = {'CA3 pyr','CA3 int','CA1 pyr','CA1 int'}

isv = 1;

foi = [3 4 5];


frCntr = 0;
for iFreq = foi
    frCntr = frCntr + 1;

    % iFreq = 3;
    % for iFreq = 1 : nfreq
    tsv = spkVec.(caseName)(:,iFreq, isv);

    for kType = 3 : nType
    
        % subplot2d(nR,nC, );
        k = frCntr;
        subplot2d(nR,nC, [sdi+k], kType);
        polarplot(angle(tsv(Sig.cell_label == kType)), abs(tsv(Sig.cell_label == kType)), ...
                  '.', 'color',cellCol(kType,:), 'markersize',10);
        % hold on
        if iFreq == foi(1), title(cellLabels{kType}); end
    end
    % legend('CA3 pyr','CA3 int','CA1 pyr','CA1 int','location','east')
    % set(gca, 'fontsize', 12)
end

%% neuron multi compartment
subplot2d(nR,nC, [sdi+1 sdi+2], 1)
% subplot2d(nR,nC, [sdi+2 sdi+3], 3)
box on

%% LFP vec
subplot2d(nR,nC, [sdi+1 sdi+2], 2)

frCntr = 0;
for iFreq = foi
    frCntr = frCntr + 1;
    tsv = lfpVec.(caseName)(:,iFreq, isv);
    plot(real(tsv), (1:32), 'linewidth',2)
    hold on
end

legend(freqBandsLab(foi), 'location','southwest')

axis tight

ax = gca;
ax.YDir = 'reverse'
