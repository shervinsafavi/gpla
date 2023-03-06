% main difference with v3 is we have stats
% the visual aestics

%%
load(fullfile(dirPath_gpla.data, 'handies', 'hpcSmltHandyVars_wAnalyticLfp.mat'))
load(fullfile(dirPath_gpla.data, 'dataHPsim', 'Session(full)_10_CA3CA1Network.mat'))

sc_threshold = 500;
nfreq = numel(freqBands);

%%
sameElecCheckInfo = [];

% statTestInfo = []; % no stat test will be applied
% new stat parameters
statTestInfo.testType = 'spike-jittering';
statTestInfo.SVspectrumStatsType = 'RMT-heuristic';
statTestInfo.jitterType = 'group-preserved-interval-jittering';

statTestInfo.nJtr = 500;
statTestInfo.alphaValue = 0.05;
statTestInfo.spkSF = Fs;


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

% commented out as it the results is saved 

% tic
% for iFreq = 1 : nfreq
%     lfpPh = lfpAnalSig_allSepFreq{iFreq};
%    
%     statTestInfo.jitterWinWidth = 1 / mean(freqBands{iFreq});
% 
%     [lfpVec.(caseName)(:,iFreq), spkVec.(caseName)(:,iFreq), gPLV.(caseName)(iFreq), ...
%      ~, allStats.(caseName)(iFreq), rawSvdStuff.(caseName)(iFreq)] = ...
%         tngpla(spkTrains, lfpPh ,[], sc_threshold, regionSpkInd, [] , 1, ...
%                statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);
% end
% 
% % collect stuff for significance
% for ifreq = 1 : nfreq
%    
%     empSig(ifreq) = allStats.(caseName)(ifreq).gPLV_stats.nullHypoReject;
    
%     dimRatio = size(rawSvdStuff.(caseName)(ifreq).couplingMatrix, 1) ...
%         / size(rawSvdStuff.(caseName)(ifreq).couplingMatrix, 2);
%     lambda(ifreq) = (1 + dimRatio^.5) ^ 2;
%     sgplv(ifreq) = gPLV.(caseName)(ifreq) ^ 2 / ...
%         size(rawSvdStuff.(caseName)(ifreq).couplingMatrix, 2);
% 
% end
% 
% theoSig = sgplv > lambda;


% % save the computed data
% save(fullfile('dat', 'ca1_gpla_v4_partialVars.mat'), ...
%      'empSig', 'theoSig', 'lfpVec', 'spkVec', 'gPLV', ... 
%      'allStats', 'rawSvdStuff', 'statTestInfo', ...
%      '-v7.3')


load(fullfile('dat', 'ca1_gpla_v4_partialVars.mat')); 

freqBandsLab = {...
    '1-5';...
    '5-20'; ...
    '20-40'; ...
    '40-80'; ...
    '80-180' ...
    };


%% plot gPLV
% subplot2d(nR,nC, sdi+1, [1 2]);

% delete(sbh)

k = 3;
[~, sbh] = subplot2d(nR,nC, k, [3 4]-2);
relocateFigureStuff(sbh, [pos_lastRowLfpDemo(1) NaN NaN NaN])
% displaceFigureStuff(sbh, [NaN .025 NaN NaN])
% displaceFigureStuff(sbh, [NaN .045 NaN NaN])
displaceFigureStuff(sbh, [.04 .045 NaN NaN])

% bar(gPLV.(caseName), 'k')

mo = 30; % vertical offset for the trianglesl
vc.f11.sigTriangleMS = 5;
% vc.f11.sigTriangleLW = 1.1;
vc.f11.sigTriangleLW = 1.5;
vc.f11.sigTriangleLFS = 5.5; % legend font size
vc.f11.lfs = 6;

vc.f11.gPLVxTickFS = 6.5;

hold all

triagnelXos = .15; % triagnel x offset arount the bar center
for ifreq = 1 : nFreq
    tv = gPLV.(caseName)(ifreq);
    if empSig(ifreq) == 1, 
        plot(ifreq-triagnelXos, mo+ tv, 'v', 'color', [.7 0 .7], ...
                                'markersize',vc.f11.sigTriangleMS, ...
                                'linewidth',vc.f11.sigTriangleLW);
    end
    if  theoSig(ifreq) == 1, 
        plot(ifreq+triagnelXos, mo+ tv, 'v', 'color', [0 .7 0], ...
             'markersize',vc.f11.sigTriangleMS, ...
             'linewidth',vc.f11.sigTriangleLW);
    end
end

bar(gPLV.(caseName),'k', 'barwidth', .6);
% xlim([0 nFreq] + .5)
set(gca,'xtick', (1:nFreq))
set(gca,'xticklabel',freqBandsLab)

ylim([0 200])
set(gca, 'ytick', [0 100 200]);
set(gca, 'TickDir', 'out');
% set(gca, 'TickLength', [0.0155 0.045]);
set(gca, 'TickLength', [0.03 0.035]);
xtickangle(45)

xtickangle(45)
xlabel('Frequency band [Hz]')
ylabel('gPLV')
set(gca, 'fontSize', vc.f11.gPLVxTickFS)

[lh icons] = legend('Emprical test', 'Theoretical test', ...
            'location','northwest');

set(lh,...
    'box', 'off', ...
    'fontsize', vc.f11.sigTriangleLFS);
displaceFigureStuff(lh, [-.05 .005 NaN NaN]);
bof = +.1; % backward offset
for k = 3 : 6
    tmpOrigPos = icons(k).XData;
    icons(k).XData = tmpOrigPos + bof; 
    icons(k).MarkerSize = 5.5;
    icons(k).LineWidth = 1.5;
end


for k = 1 : 2
    icons(k).FontSize = vc.f11.lfs;
end


%% ~ plot spike vectors
cell2mat(freqBands')

nType = 4;

cellLabels = {'CA3 pyr','CA3 int','CA1 pyr','CA1 int'}

isv = 1;

foi = [3 4 5];

% off set of freq for spike vector polrar plots
ros(1) = .11;
ros(2) = .3;
ros(3) = ros(2);



frCntr = 0;
for iFreq = foi
    frCntr = frCntr + 1;

    % iFreq = 3;
    % for iFreq = 1 : nfreq
    tsv = spkVec.(caseName)(:,iFreq, isv);

    for kType = 3 : nType
    
        % subplot2d(nR,nC, );
        k = frCntr;

        [~, sbh] = subplot2d(nR,nC, [sdi+k], kType);
        % displaceFigureStuff(sbh, [NaN NaN .005 .005])
        displaceFigureStuff(sbh, [NaN -.0295 .005 .005]);

        if kType == 4
            displaceFigureStuff(sbh, [-.035 NaN NaN NaN]);
        end
        
        polarplot(angle(tsv(Sig.cell_label == kType)), abs(tsv(Sig.cell_label == kType)), ...
                  '.', 'color',cellCol(kType,:), 'markersize',8);
        bfypp1

        set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
        set(gca, 'fontsize', vc.f11.gPLVxTickFS-1)
        % hold on
        if iFreq == foi(1), 
            th = title(cellLabels{kType}); 
            set(th, 'fontweight', 'normal')
            set(th, 'fontsize', 9)
        end
    end

    if kType == 4
        rRange = get(gca, 'rlim')
        text(0, rRange(2) + ros(frCntr), [freqBandsLab{iFreq}, 'Hz'], ...
             'fontsize', 8.5)
    end

    % legend('CA3 pyr','CA3 int','CA1 pyr','CA1 int','location','east')
    % set(gca, 'fontsize', 5)
end

%% neuron multi compartment
teh = 0.02; % temp enlarge height
[~, sbh] = subplot2d(nR,nC, [sdi+1 sdi+2], 1);
displaceFigureStuff(sbh, [NaN NaN NaN teh]);
placeHolder
% subplot2d(nR,nC, [sdi+2 sdi+3], 3)
% box on

%% LFP vec
[~, sbh] = subplot2d(nR,nC, [sdi+1 sdi+2], 2);
displaceFigureStuff(sbh, [NaN NaN -.06 teh])
frCntr = 0;
for iFreq = foi
    frCntr = frCntr + 1;
    tsv = lfpVec.(caseName)(:,iFreq, isv);
    plot(real(tsv), (1:32), 'linewidth',2)
    hold on
end

ylabel('Channel ID')
xlabel('Magnitude')
set(gca, 'fontSize', vc.f11.gPLVxTickFS)

% [lh, icons] = legend(freqBandsLab(foi), 'location','southwest');
% set(lh, 'box', 'off')
% displaceFigureStuff(lh, [NaN .045 0 0])

[lh, icons] = legend(freqBandsLab(foi), 'location','northwest');
set(lh, 'box', 'off')
displaceFigureStuff(lh, [.03 -.01 0 0])

axis tight
lvl = .5; % LFP vector lim
xlim(lvl * [-1 1])
% xlim([-.5 0])

% set(gca, 'TickDir', 'out');
set(gca, ...
    'TickLength', [0.03 0.035], ...
    'TickDir',    'out', ...
    'ytick',      (8:8:32), ...
    'xtick',      [-lvl 0 lvl] ...
    );

% set(gca, 'TickDir', 'out');
% set(gca, 'TickLength', [0.0155 0.045]);
% set(gca, 'TickLength', [0.03 0.035]);

% grid on

bof = 0.3; % backward offset
for k = 4 : 2: 9
    tmpOrigPos = icons(k).XData(2);
    icons(k).XData(2) = tmpOrigPos - bof; 
end

for k = 1 : 3
    icons(k).Position(1) = 0.3;
end

ax = gca;
ax.YDir = 'reverse'
% ax.XDir = 'reverse'
box off