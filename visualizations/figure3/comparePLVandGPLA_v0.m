% this is supposed to to only the computtiaon, but it also contain
% the plot stufff

%% 
% this snipt stolen from v0.m from line around 219 

% clear all
% ignit_gpla
% 
pds = ignit();
run('../vizConventions')

%%
caseType = 'gplaBenefit';
iSV = 1; % in this simulation we only focus on the fist/largest SV 

%%
% We simulate Poisson spike trains transiently phase-locked  to an
% oscilation in a prticular freqency band.

% All simulations share the following common parameters:

signalParams = struct ...
    ( ...
        'nCh',      1, ...          % number of channel with LFP
        'nUnit',    18, ...          % number of spiking units
        'SF',       1e3, ... % Hz    % sampling frequency
        'nTr',      1, ...          % number of trials
        'signalLength', ...  % sec   % duration of the signals
        10 ...
        );

spikeTrainParams = struct ...
    (...
        'avefiringRate', ...  % Hz   % average firing rate in the transient epochs
        20 ...
        );

globalDynamicsParams = struct ...
    (...
        'oscFreq',  22.5, ... % Hz   % oscillation frequency for LFP and phase-locking
        'nCycl',    20, ...          % number oscillation cycles each transient has
        'syncSigProportion', ...     % the proportion/fraction of the signal with transients
        .8, ...   
        'lfpPhaseNoise_kappa', ...   % concentration parameter $\kappa$ for noise on LFP phase
        10, ...
        'whiteNoise_sigma', ...      % variance for white noise for non-oscillatory periods
        0.05 ...  
        );

% refPhase = pi/2;
% phaseShifter = [refPhase; refPhase+(1*2*pi/3); refPhase+(2*2*pi/3)];

% couplingParams = struct ...
%     (...
%     'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
%                 0.4 , ...      % this specific value will ensure that theoritical phase locking value will be 0.3   
%     'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
%          phaseShifter(:) ... % rad  
%     );

% whiteNoise_sigma = 0.05; 

%% Transients epochs
% Transients epochs are specified by Poisson events. The shoule have
% certain IEI 

transientUnit_duration = (globalDynamicsParams.nCycl * (1/globalDynamicsParams.oscFreq));
for iTr = 1 : signalParams.nTr
    eventCount = round(globalDynamicsParams.syncSigProportion * signalParams.signalLength) ...
        / transientUnit_duration;
    eventRate = eventCount / signalParams.signalLength;
    iEvent = 1;
    eventTrain{iTr}(iEvent) = -1 * log(rand(1)) / eventRate;
    while eventTrain{iTr}(iEvent) < signalParams.signalLength
        tmp_iei = (-1 * log(rand(1)) / eventRate);
        if (tmp_iei >= transientUnit_duration)
            iEvent = iEvent + 1;
            eventTrain{iTr}(iEvent) = eventTrain{iTr}(iEvent - 1) + tmp_iei;
        end
    end
    
    eventTrain{iTr} = eventTrain{iTr}(eventTrain{iTr} < signalParams.signalLength - transientUnit_duration);
end

%% Filter 
% freqBand = [15 30];
freqBand = globalDynamicsParams.oscFreq + [-2.5 +2.5];
wn = freqBand / (signalParams.SF/2);
filterOrder = 2;
[b, a] = butter(filterOrder, wn, 'bandpass');

%% define stat info
statTestInfo.testType = 'spike-jittering';
statTestInfo.SVspectrumStatsType = 'RMT-heuristic';
statTestInfo.jitterType = 'interval-jittering';

statTestInfo.nJtr = 100;
statTestInfo.alphaValue = 0.05;
statTestInfo.spkSF = signalParams.SF;

statTestInfo.jitterWinWidth = 1 / mean(globalDynamicsParams.oscFreq);

sameElecCheckInfo = [];

%% viz conve
iTr_forViz = 1;
iEvent_forViz = 2;
tmpOffSet = 100;
transDur = transientUnit_duration * signalParams.SF;


%%
nCase = 1;


%% Case III: fig 1C
% Similar to case I, but theres is a phase shift between successive
% channels/units.

% define parameters
iCase = 1;
nChunk = 3;
% phaseShifter = zeros(nChunk, signalParams.nUnit / nChunk);
pPhShif = zeros(signalParams.nUnit / nChunk, nChunk); % pre phase shifter
                                                      % phaseShifter = zeros(signalParams.nUnit, 1);
                                                      % chankEdge = linspace(1, signalParams.nUnit, nChank);
                                                      % tmpInd = (reshape((1:signalParams.nUnit), signalParams.nUnit / nChunk, nChunk))';
for iChunk = 1 : nChunk
    pPhShif(1 : signalParams.nUnit / nChunk, iChunk) = iChunk *  2*pi/nChunk;
end
phaseShifter = pPhShif(:, randperm(nChunk));
% phaseShifter = linspace(0, pi, signalParams.nCh)';

couplingParams(iCase) = struct ...
    (...
        'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
        5 , ...
        'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
        phaseShifter(:) ... % rad  
        );

% .3 is a value which sometimes lead to significant locling for PLV
% sometimes not

% run the simulation 
[lfpLikeSig.(caseType){iCase}, spikeTrains.(caseType){iCase}] = ...
    smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, ...
                                couplingParams(iCase), signalParams, eventTrain);

for iTr = 1 : signalParams.nTr
    for iCh = 1 : signalParams.nCh
        tmpLfp = lfpLikeSig.(caseType){iCase}(iCh,:, iTr);
        tmpLfp_filtered = filtfilt(b, a, tmpLfp);
        tmpAnlcLfp = hilbert(tmpLfp_filtered);
        lfpLikePhase.(caseType){iCase}(iCh,:, iTr) = angle(tmpAnlcLfp);
        lfpLikeAnalSig.(caseType){iCase}(iCh,:, iTr) = tmpAnlcLfp;
    end
end



% [lfpVec.mthdIllus(:,iCase), spkVec.mthdIllus(:,iCase), gPLV.mthdIllus(iCase), ...
%  ~, gPLV_stats.mthdIllus(iCase), PLV.mthdIllus(:,:, iCase), ~, PLV_stats.mthdIllus(iCase)] = ...
%     gpla(spikeTrains.mthdIllus{iCase}, lfpLikePhase.mthdIllus{iCase}, [],[],[],[],[], statTestInfo);

% [lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
%  ~, ~, rawSvdStuff.(caseType)(iCase)] = ...
%     tngpla(spikeTrains.(caseType){iCase}, lfpLikeAnalSig.(caseType){iCase}, [], [], [], [] , [], ...
%            [], iSV, sameElecCheckInfo, 'nSpk-square-root', 2);

% [lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
%  ~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
%     tngpla(spikeTrains.(caseType){iCase}, lfpLikePhase.(caseType){iCase}, [], [], [], [] , [], ...
%            statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 0);


%%
iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpRng = tmpInd-tmpOffSet : tmpInd+transDur+tmpOffSet ;

%%
nR = 4; nC = 4;

%% Publication figure
% clf
iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpOffSet = 70;

tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpMid = round(eventTrain{iTr}(iEvent) * signalParams.SF ...
               + transDur / 2);
tmpRng = tmpMid-tmpOffSet : tmpMid+tmpOffSet ;


%%

% subplot2(nR,nC, 1,1)
% plot(lfpLikeSig.(caseType){iCase}(iCh,tmpInd-100:tmpInd+4000 , iTr), 'color','k'); %colorbar
% axis tight
% ylim([-1.5 1.5])

%%
% subplot2(nR,nC, 1, 2)

% subplot2(nR,nC, 2, 1)
% plot_EveConData(lfpLikeSig.(caseType){iCase}(1,tmpRng, iTr), full(spikeTrains.(caseType){iCase}{iTr}(:,tmpRng)), ...
%                 {'color','k', 'linewidth',vc.lw_lfp}, {'color','r', 'linewidth',vc.lw_spk})        
% % for bar plot 
% % pplt.CData(iCase,:) = iCase*(0.5/nCase)*vc.illusCases;

%%
% subplot2(nR,nC, 2,2);
% polarplot(angle(spkVec.(caseType)(:,iCase)), abs(spkVec.(caseType)(:,iCase)), ...
%           'color',vc.spkv, 'LineStyle','none', 'marker','.')
% rlim([0 .5])


%% compare PLV and GPLA
% new parameters 
couplingParams(iCase).lockingStrength_kappa = .3;
signalParams.nUnit = 3;
tv = linspace(0, 2*pi, signalParams.nUnit+1);
phaseShifter = tv(1:end-1);

couplingParams(iCase) = struct ...
    (...
        'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
        .3, ...
        'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
        phaseShifter(:) ... % rad  
        );


% for iChunk = 1 : nChunk
%     pPhShif(1 : signalParams.nUnit / nChunk, iChunk) = iChunk *  2*pi/nChunk;
% end
% phaseShifter = pPhShif(:, randperm(nChunk));

% nRel = 500;
nRel = 5000;

clear lfpLikePhase lfpLikeAnalSig lfpVec spkVec statSummary

for iRel = 1 : nRel
    
    % run the simulation 
    [lfpLikeSig.(caseType){iCase}, spikeTrains.(caseType){iCase}] = ...
        smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, ...
                                    couplingParams(iCase), signalParams, eventTrain);

    for iTr = 1 : signalParams.nTr
        for iCh = 1 : signalParams.nCh
            tmpLfp = lfpLikeSig.(caseType){iCase}(iCh,:, iTr);
            tmpLfp_filtered = filtfilt(b, a, tmpLfp);
            tmpAnlcLfp = hilbert(tmpLfp_filtered);
            lfpLikePhase.(caseType){iCase}(iCh,:, iTr) = angle(tmpAnlcLfp);
            lfpLikeAnalSig.(caseType){iCase}(iCh,:, iTr) = tmpAnlcLfp;
        end
    end

    % for gPLV
    [lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
     ~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
        tngpla(spikeTrains.(caseType){iCase}, lfpLikeAnalSig.(caseType){iCase}, [], [], [], [] , [], ...
               statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', ...
               2);
    
    statSummary.gPLV(iRel) = allStats.gplaBenefit.gPLV_stats.nullHypoReject;
    
    % for PLV
    % we use GPLA code, but without whitening and use usual normalization
    [~,~,~,~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
        tngpla(spikeTrains.(caseType){iCase}, lfpLikePhase.(caseType){iCase}, [], [], [], [] , [], ...
               statTestInfo, iSV, sameElecCheckInfo, 'nSpk', ...
               0);
    
    PLVexamples(iRel, :) = rawSvdStuff.(caseType)(iCase).couplingMatrix;
        statSummary.PLV(iRel, :) = allStats.gplaBenefit.PLV_stats.nullHypoReject;
    
    % for pPLV
    [pPLV.(caseType)(iCase), ~, ~, pPLV_stat.(caseType)(iCase)] = ...
    ppla(spikeTrains.(caseType){iCase}, lfpLikePhase.(caseType){iCase}, [],[],[], statTestInfo);
    statSummary.pPLV(iRel) = pPLV_stat.(caseType)(iCase).nullHypoReject;

end

%%
barData(1 : signalParams.nUnit) = 100 * sum(statSummary.PLV) / nRel;
barData(signalParams.nUnit + 1) = 100 * sum(statSummary.pPLV) / nRel;
barData(signalParams.nUnit + 2) = 100 * sum(statSummary.gPLV) / nRel;


%%
[~,~,~,~,~, rawSvdStuff.(caseType)(iCase)] = ...
        tngpla(spikeTrains.(caseType){iCase}, lfpLikePhase.(caseType){iCase}, [], [], [], [] , [], ...
               [], iSV, sameElecCheckInfo, 'nSpk', 0);

subplot2(nR,nC, 1, 3);

for iRel = 1 : nRel
    tmp = PLVexamples(iRel, :);
    polarplot(angle(tmp), abs(tmp),...
              '.', 'color',.6*ones(1,3), 'MarkerSize', vc.f7.bulletSize/5); 
    hold on
end


% tmp = rawSvdStuff.gplaBenefit.couplingMatrix;
tmp = mean(PLVexamples, 1);
polarplot(angle(tmp), abs(tmp),...
    'k.', 'MarkerSize', vc.f7.bulletSize);

rlim([0 .5])
ax = gca;
ax.LineWidth = 1.5;
rruler = ax.RAxis;
rruler.Label.String = 'PLV';
ax.RAxisLocation = 45;

%%
subplot2(nR,nC, 1, 4);
bar(barData, vc.f7.barWidth, 'facecolor', 'k')
grid on

ylabel('Successful detection [%]')

% xticklabels({'Unit 1 PLV','Unit 2 PLV','Unit 3 PLV', 'pPLV', 'gPLV'}); %
% xticklabels({'Unit 1','Unit 2','Unit 3', 'pPLV', 'gPLV'});
xticklabels({'U 1','U 2','U 3', 'pPLV', 'gPLV'});
% xticklabels({'','PLV','Unit 3', 'pPLV', 'gPLV'});
xtickangle(45)

%% save 
% we save the entire workspace 
save(fullfile('dat','comparePLVandGPLA_v0'));
