% clear all
% ignit_gpla

vc.f5.gPLVxTickFS = 4;

%% load data with analytical signals and whiten them
load(fullfile(dirPath_gpla.data, 'handies', ...
              'uate2010handyVars_wAnalyticLfp.mat'))

nfreq = size(freqBands,1);

sameElecCheckInfo = [];
iSV = 1;

caseName = 'uate2010'; 


% commented out as it the results is saved 
% % new stat parameters
% statTestInfo.testType = 'spike-jittering';
% statTestInfo.SVspectrumStatsType = 'RMT-heuristic';
% statTestInfo.jitterType = 'group-preserved-interval-jittering';

% statTestInfo.nJtr = 100;
% statTestInfo.alphaValue = 0.05;
% statTestInfo.spkSF = Fs;

% %% GPLA


% for iFreq = 1 : nfreq
%     freqBand = freqBands(iFreq,:);
%     lfpPh = lfpAnalSig_allSepFreq{iFreq};  

%     statTestInfo.jitterWinWidth = 1 / mean(freqBands(iFreq,:));

%     [lfpVec.(caseName){iFreq}, spkVec.(caseName){iFreq}, gPLV.(caseName)(iFreq), ...
%      ~, allStats.(caseName)(iFreq), rawSvdStuff.(caseName)(iFreq)] = ...
%         tngpla(sparseSPT, lfpPh ,[], sc_threshold, [], [] , 1, ...
%                statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);
% end

% %% collect stuff for significance
% for ifreq = 1 : nfreq
    
%     empSig(ifreq) = allStats.(caseName)(ifreq).gPLV_stats.nullHypoReject;
    
%     dimRatio = size(rawSvdStuff.(caseName)(ifreq).couplingMatrix, 1) / size(rawSvdStuff.(caseName)(ifreq).couplingMatrix, 2);
%     lambda(ifreq) = (1 + dimRatio^.5) ^ 2;
%     sgplv(ifreq) = gPLV.(caseName)(ifreq) ^ 2 / ...
%         size(rawSvdStuff.(caseName)(ifreq).couplingMatrix, 2);
    
% end

% theoSig = sgplv > lambda

% instead load 
load(fullfile('dat', 'gplaAcrossFreqs_v0.mat')); 

%% save the computed data
% save(fullfile('dat', 'gplaAcrossFreqs_v0.mat'), ...
%      'empSig', 'theoSig', 'lfpVec', 'spkVec', 'gPLV', ... 
%      'allStats', 'rawSvdStuff', 'statTestInfo', ...
%      '-v7.3')

%% plot
hold all
for ifreq = 1 : 8
    tv = gPLV.(caseName)(ifreq);
    if empSig(ifreq) == 1, plot(ifreq-.2, 14+ tv, 'bv'); end
    if theoSig(ifreq) == 1, plot(ifreq+.2, 14+ tv, 'rv'), end
end


bar(gPLV.(caseName),'k');
set(gca,'xtick', (1:8))
set(gca,'xticklabel',freqBandsLab)
xtickangle(45)
xlabel('Frequency band (Hz)')
ylabel('gPLV')
set(gca, 'fontSize', 14)
grid on

legend('Emprical test', 'Theoretical test')

%% 