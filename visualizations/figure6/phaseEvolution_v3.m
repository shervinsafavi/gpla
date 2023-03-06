% difference with v2 is in the plotting

% we plot the fine grained phase evo

load(fullfile(pds.src, 'explorations', '261_shp_fineGrainedHighFreq', 'dat', ...
     'hpcSmltHandyVars_wAnalyticLfp.mat'))

nfreq = numel(lfpAnalSig_allSepFreq);

lfpAnalSig_allSepFreq_old = lfpAnalSig_allSepFreq;
clear lfpAnalSig_allSepFreq
for ifreq = 1 : nfreq
    lfpAnalSig_allSepFreq{ifreq} = lfpAnalSig_allSepFreq_old{ifreq}(regionLfpsInd, :,:);
end
n.nCh = 32;

clear lfpAnalSig_allSepFreq_old

%% GPLA
statTestInfo = [];
caseName = 'ca1anal_fineGrainedHighFreq';
% clear lfpVec spkVec
iSV = 1; % this will change in later part of the analysis

for iFreq = 1 : nfreq
    lfpPh = lfpAnalSig_allSepFreq{iFreq};

    [lfpVec.(caseName)(:,iFreq), spkVec.(caseName)(:,iFreq), gPLV.(caseName)(iFreq), ...
     ~, ~, rawSvdStuff.(caseName)(iFreq)] = ...
        tngpla(spkTrains, lfpPh ,[], sc_threshold, regionSpkInd, [] , 1, ...
               statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);
end


%% cmpt phase spike-LFP phase difference

for iFreq = 1 : nfreq
    tlv = lfpVec.(caseName)(:,iFreq, isv);
    
    % pyr
    kType = 3;
    cioi = find(Sig.cell_label == kType);
    tsv = spkVec.(caseName)(cioi,iFreq, isv);
    % spkLfpRatio.pyr(iFreq) = tsv / tlv(16);
    % spkLfpRatio(iFreq, 1) = angle(nanmean(tsv / tlv(16)));
    % spkLfpRatio(iFreq, 2) = unwrap(pi + angle(nanmean(tsv / tlv(16))));
    % spkLfpRatio(iFreq, 2) = (pi + angle(nanmean(tsv / tlv(16))));
    spkLfpRatio(iFreq, 1) = wrapTo2Pi(angle(nanmean(tsv / tlv(16))));
    
    % int
    kType = 4;
    cioi = find(Sig.cell_label == kType);
    tsv = spkVec.(caseName)(cioi,iFreq, isv);
    % spkLfpRatio.int(iFreq) = tsv / tlv(16);
    % spkLfpRatio(iFreq, 2) = angle(nanmean(tsv / tlv(16)));
    % spkLfpRatio(iFreq, 2) = unwrap(pi + angle(nanmean(tsv / tlv(16))));
    % spkLfpRatio(iFreq, 2) = (pi + angle(nanmean(tsv / tlv(16))));
    spkLfpRatio(iFreq, 2) = wrapTo2Pi(angle(nanmean(tsv / tlv(16))));
end

%%  plot
% fh = bar(spkLfpRatio,'FaceColor','flat');

freqCenters = mean(cell2mat(freqBands'), 2);


% fh = bar(spkLfpRatio,'FaceColor','flat');

% for k = 1:size(spkLfpRatio,2)
%     semilogx(freqCenters, spkLfpRatio(:, k), 'linewidth', 2, ...
%              'color',cellCol(k+2, :))
%     hold on
%     % fh(k).CData = cellCol(k+2, :);
% end

for k = 1:size(spkLfpRatio,2)
    semilogx(freqCenters, rad2deg(spkLfpRatio(:, k)), 'linewidth', 2, ...
             'color',cellCol(k+2, :))
    hold on
    % fh(k).CData = cellCol(k+2, :);
end

box off
xlabel('Frequency [Hz]')
ylabel('\Delta\phi [deg]')
ylim([0 360])

%
set(gca, ...
    'TickLength', [0.03 0.035], ...
    'TickDir',    'out', ...
    'ytick',      [0 180 360] ...
    );


%%
lh = legend(cellLabels{3}, cellLabels{4}, 'location','southwest');
set(lh, 'box', 'off')
displaceFigureStuff(lh, [-.01 -.005 0 0])



% ax = gca;
% ax.YDir = 'reverse'

