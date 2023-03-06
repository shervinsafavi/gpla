% different in font size comapre to v0

% clf
% [~, sbh] = subplot2(nR,nC, nred+1, [1 2]);

% % 
% clear all
% ignit_gpla
pds = ignit()

vc.f5.gPLVxTickFS = 5;
mo = 26; % vertical offset for the trianglesl

% % load data with analytical signals and whiten them
% *** check if need to be uncommeted
% load(fullfile(dirPath_gpla.data, 'handies', ...
%               'uate2010handyVars_wAnalyticLfp.mat'))

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

% % save the computed data
% save(fullfile('dat', 'gplaAcrossFreqs_v0.mat'), ...
%      'empSig', 'theoSig', 'lfpVec', 'spkVec', 'gPLV', ... 
%      'allStats', 'rawSvdStuff', 'statTestInfo', ...
%      '-v7.3')

% % plot

% clf
% subplot2(nR,nC, nred+1, [1 2]);

mo = 44; % vertical offset for the trianglesl
vc.f5.sigTriangleMS = 4;
vc.f5.sigTriangleLW = 1.1;
vc.f5.sigTriangleLFS = 5.5; % legend font size

hold all
for ifreq = 1 : 8
    tv = gPLV.(caseName)(ifreq);
    if empSig(ifreq) == 1, plot(ifreq-.2, mo+ tv, 'v', 'color', [.7 0 .7], ...
                                'markersize',vc.f5.sigTriangleMS, ...
                                'linewidth',vc.f5.sigTriangleLW);
    end
    if  theoSig(ifreq) == 1, 
        plot(ifreq+.2, mo+ tv, 'v', 'color', [0 .7 0], ...
             'markersize',vc.f5.sigTriangleMS, ...
             'linewidth',vc.f5.sigTriangleLW);
    end
end


% pos = get(sbh, 'Position')
% posx = pos(1);
% posy = pos(2) + 0.03;
% set(sbh,...
%     'Position', [posx posy pos(3) pos(4)])


bar(gPLV.(caseName),'k');
set(gca,'xtick', (1:8))
set(gca,'xticklabel',freqBandsLab)

%
ylim([0 300])
set(gca, 'ytick', [0 150 300]);
set(gca, 'TickDir', 'out');
% set(gca, 'TickLength', [0.0155 0.045]);
set(gca, 'TickLength', [0.0355 0.035]);
xtickangle(45)
% ax.TickLength = [0.0155 0.035];

%%

xtickangle(45)
xlabel('Frequency band [Hz]')
ylabel('gPLV')
set(gca, 'fontSize', vc.f5.gPLVxTickFS)
% grid on

% lh = legend('Emprical test', 'Theoretical test', ...
%             'location','northeast');
% displaceFigureStuff(lh, [0.05 .041 0 0]);
% set(lh,...
%     'box', 'on', ...
%     'fontsize', vc.f5.sigTriangleLFS);

% LEGEND ON TOP RIGHT CORNER
% lh = legend('Emprical test', 'Theoretical test', ...
%             'location','northeast');
% displaceFigureStuff(lh, [0.05 .041 0 0]);
% set(lh,...
%     'box', 'on', ...
%     'fontsize', vc.f5.sigTriangleLFS);
[lh] = legend('Emprical test', 'Theoretical test', ...
            'location','northwest');
displaceFigureStuff(lh, [-0.02 .061 0 0]);
set(lh,...
    'box', 'on', ...
    'fontsize', vc.f5.sigTriangleLFS);

% for k = 1 : 2
%     icons(k).FontSize = vc.f5.lfs;
% end


% h = legend('Emprical test', 'Theoretical test')
% pos = get(h,'Position')

% % top right
% posx = 0.3056;
% posy = 0.4640;
% legWidth = 0;

% set(h,...
%     'Position', [posx posy legWidth pos(4)], ...
%     'box', 'on', ...
%     'fontsize', vc.f5.sigTriangleLFS);

