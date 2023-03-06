%% Figure 2: Illustrative simulations
% This figure depict the application of the method on simple simulations

%% initiatation 
% (e.g add necessary packages and functions)
clear all
clf
ignit;

%% prepare the figure

% visualization conventions
vc = get_vizConventions();


fh = gcf;

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f1.w;
fh.Position(4)  = vc.f1.h;

% number of rows and columns in the subplot
nR = 5;
nC = 4;


%% simulation parameters 
caseType = 'mthdIllus';
iSV = 1; % in this simulation we only focus on the fist/largest SV

%%
% We simulate Poisson spike trains transiently phase-locked (or not) to an
% oscillation in a particular freq ency band.

% All simulations share the following common parameters:
signalParams = struct ...
    ( ...
        'nCh',      1, ...          % number of channel with LFP
        'nUnit',    18, ...          % number of spiking units
        'SF',       1e3, ... % Hz    % sampling frequency
        'nTr',      10, ...          % number of trials
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

% whiteNoise_sigma = 0.05; 

%%
nCase = 4;

%% Transients epochs
% Transients epochs are specified by Poisson events. The shoule have
% certain inter-event-interval (IEI)

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

%% define statistics-testing information
statTestInfo.testType = 'spike-jittering';
statTestInfo.SVspectrumStatsType = 'RMT-heuristic';
statTestInfo.jitterType = 'group-preserved-interval-jittering';

statTestInfo.nJtr = 1;
statTestInfo.alphaValue = 0.05;
statTestInfo.spkSF = signalParams.SF;

statTestInfo.jitterWinWidth = 1 / mean(globalDynamicsParams.oscFreq);

sameElecCheckInfo = [];

%% some conventions for visualization
iTr_forViz = 1;
iEvent_forViz = 2;
tmpOffSet = 100;
transDur = transientUnit_duration * signalParams.SF;

%% Case I
% In this case, neurons are considered to be homogeneous, i.e.
% modulation is identical across all neurons.

% define parameters
iCase = 1;
couplingParams(iCase) = struct ...
    (...
        'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
        10 , ...
        'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
        pi ...% rad  
        );

% run the simulation 
[lfpLikeSig.mthdIllus{iCase}, spikeTrains.mthdIllus{iCase}] = ...
    smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, ...
                                couplingParams(iCase), signalParams, eventTrain);

for iTr = 1 : signalParams.nTr
    for iCh = 1 : signalParams.nCh
        tmpLfp = lfpLikeSig.mthdIllus{iCase}(iCh,:, iTr);
        tmpLfp_filtered = filtfilt(b, a, tmpLfp);
        tmpAnlcLfp = hilbert(tmpLfp_filtered);
        lfpLikePhase.mthdIllus{iCase}(iCh,:, iTr) = angle(tmpAnlcLfp);
        lfpLikeAnalSig.mthdIllus{iCase}(iCh,:, iTr) = tmpAnlcLfp;
    end
end

[lfpLikeSig.mthdIllus{iCase}, spikeTrains.mthdIllus{iCase}] = ...
    smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, ...
                                couplingParams(iCase), signalParams, eventTrain);


[lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
 ~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
    tngpla(spikeTrains.mthdIllus{iCase}, lfpLikeAnalSig.mthdIllus{iCase}, [], [], [], [] , [], ...
           statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);


iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpRng = tmpInd-tmpOffSet : tmpInd+transDur+tmpOffSet ;

%% Case II
% Similar to case I, but there is a phase shift between successive
% channels/units.

% define parameters
iCase = 2;
phaseShifter = linspace(0, pi, signalParams.nUnit)';

couplingParams(iCase) = struct ...
    (...
        'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
        10 , ...
        'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
        phaseShifter ... % rad  
        );

% run the simulation 
[lfpLikeSig.mthdIllus{iCase}, spikeTrains.mthdIllus{iCase}] = ...
    smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, ...
                                couplingParams(iCase), signalParams, eventTrain);

for iTr = 1 : signalParams.nTr
    for iCh = 1 : signalParams.nCh
        tmpLfp = lfpLikeSig.mthdIllus{iCase}(iCh,:, iTr);
        tmpLfp_filtered = filtfilt(b, a, tmpLfp);
        tmpAnlcLfp = hilbert(tmpLfp_filtered);
        lfpLikePhase.mthdIllus{iCase}(iCh,:, iTr) = angle(tmpAnlcLfp);
        lfpLikeAnalSig.mthdIllus{iCase}(iCh,:, iTr) = tmpAnlcLfp;

    end
end

[lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
 ~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
    tngpla(spikeTrains.mthdIllus{iCase}, lfpLikeAnalSig.mthdIllus{iCase}, [], [], [], [] , [], ...
           statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);

iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpRng = tmpInd-tmpOffSet : tmpInd+transDur+tmpOffSet ;


%% Case III
% Similar to case I, but theres is a phase shift between successive
% channels/units.

% define parameters
iCase = 3;
nChunk = 3;
pPhShif = zeros(signalParams.nUnit / nChunk, nChunk); % pre phase shifter
                                                      % phaseShifter = zeros(signalParams.nUnit, 1);
                                                      % chankEdge = linspace(1, signalParams.nUnit, nChank);
                                                      % tmpInd = (reshape((1:signalParams.nUnit), signalParams.nUnit / nChunk, nChunk))';
for iChunk = 1 : nChunk
    pPhShif(1 : signalParams.nUnit / nChunk, iChunk) = iChunk *  2*pi/nChunk;
end
phaseShifter = pPhShif(:, randperm(nChunk));

couplingParams(iCase) = struct ...
    (...
        'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
        10 , ...
        'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
        phaseShifter(:) ... % rad  
        );

% run the simulation 
[lfpLikeSig.mthdIllus{iCase}, spikeTrains.mthdIllus{iCase}] = ...
    smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, ...
                                couplingParams(iCase), signalParams, eventTrain);

for iTr = 1 : signalParams.nTr
    for iCh = 1 : signalParams.nCh
        tmpLfp = lfpLikeSig.mthdIllus{iCase}(iCh,:, iTr);
        tmpLfp_filtered = filtfilt(b, a, tmpLfp);
        tmpAnlcLfp = hilbert(tmpLfp_filtered);
        lfpLikePhase.mthdIllus{iCase}(iCh,:, iTr) = angle(tmpAnlcLfp);
        lfpLikeAnalSig.mthdIllus{iCase}(iCh,:, iTr) = tmpAnlcLfp;
    end
end

[lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
 ~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
    tngpla(spikeTrains.mthdIllus{iCase}, lfpLikeAnalSig.mthdIllus{iCase}, [], [], [], [] , [], ...
           statTestInfo, iSV, sameElecCheckInfo, 'nSpk-square-root', 2);

iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpRng = tmpInd-tmpOffSet : tmpInd+transDur+tmpOffSet ;

%% Case IV

% define parameters
iCase = 4;

couplingParams(iCase) = struct ...
    (...
        'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
        0, ...
        'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
        pi/2 ...  % rad  
        );

% run the simulation 
[lfpLikeSig.mthdIllus{iCase}, spikeTrains.mthdIllus{iCase}] = ...
    smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, ...
                                couplingParams(iCase), signalParams, eventTrain);

for iTr = 1 : signalParams.nTr
    for iCh = 1 : signalParams.nCh
        tmpLfp = lfpLikeSig.mthdIllus{iCase}(iCh,:, iTr);
        tmpLfp_filtered = filtfilt(b, a, tmpLfp);
        tmpAnlcLfp = hilbert(tmpLfp_filtered);
        lfpLikePhase.mthdIllus{iCase}(iCh,:, iTr) = angle(tmpAnlcLfp);
        lfpLikeAnalSig.mthdIllus{iCase}(iCh,:, iTr) = tmpAnlcLfp;
    end
end

[lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
 ~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
    tngpla(spikeTrains.mthdIllus{iCase}, lfpLikeAnalSig.mthdIllus{iCase}, [], [], [], [] , [], ...
           statTestInfo, iSV, sameElecCheckInfo, 'nSpk', 0);


iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpRng = tmpInd-tmpOffSet : tmpInd+transDur+tmpOffSet ;

%% manuscript figure
clf
iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpOffSet = 70;

tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpMid = round(eventTrain{iTr}(iEvent) * signalParams.SF ...
               + transDur / 2);
tmpRng = tmpMid-tmpOffSet : tmpMid+tmpOffSet ;

subplot2(nR,nC, 1,(1 : nC-1));

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

ttr = tmpInd-100:tmpInd+4000;
tlfp = lfpLikeSig.(caseType){iCase}(iCh,ttr , iTr);


plot(ttr, tlfp, ...
     'color','k', 'linewidth', 1); %colorbar

axis tight

tyr = 1.3 * [-1 1]; % temp Y range
ylim(tyr);


yslos = .1;

% xstos = .25;
ystos = .11;
stoslw = 1.5;
slosfs = 8.5;

plotScaleBar1(ttr, tlfp, ...
              '500 ms', NaN, 500, NaN, ...
              yslos, NaN, ystos,  NaN, ...
              stoslw, slosfs); 


axis off
box off

rectangle(gca, 'Position',[tmpRng(1) -1.1 numel(tmpRng) 2.2], ...
          'EdgeColor', 'b', 'linewidth', 1.2);


subplot2(nR,nC, 1, nC);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
pplt = bar(gPLV.mthdIllus, 'FaceColor','flat');

tym = 200;
ylim([0 tym]);
set(gca, 'ytick', [0 tym/2 tym])

ylabel('gPLV')
xticklabels({'Model 1','Model 2','Model 3', 'Model 4'});
xtickangle(45)

box off

% coarse
for iCase = 1 : nCase
    subplot2(nR,nC, iCase+1, 1);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
    axis on
end

for iCase  = 1 : nCase 
    subplot2(nR,nC, iCase+1, [2 3]);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
    plot_EveConData(lfpLikeSig.mthdIllus{iCase}(1,tmpRng, iTr), full(spikeTrains.mthdIllus{iCase}{iTr}(:,tmpRng)), ...
                    {'color',iCase*(0.5/nCase)*vc.illusCases, 'linewidth',vc.lw_lfp}, {'color','r', 'linewidth',vc.lw_spk})        
    
    % for bar plot 
    pplt.CData(iCase,:) = iCase*(0.5/nCase)*vc.illusCases;
end

%%
for iCase = 1 : nCase
    subplot2(nR,nC, iCase+1, 4);
    polarplot(angle(spkVec.mthdIllus(:,iCase)), abs(spkVec.mthdIllus(:,iCase)), ...
              'color',vc.spkv, 'LineStyle','none', 'marker','.')
    rlim([0 .5])
    
    thetaticks(0:45:315)
    tn = .5;
    rlim([0 tn])
    rticks([0 tn/2 tn]);

    ax = gca;
    ax.LineWidth = 1.5;
    rruler = ax.RAxis;
    ax.RAxisLocation = 45/2;
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
end

%%
allColors = hsv(360);
maxCoupling = 1;
nNrn = 6;

iCase = 1;
tc = allColors(180, :);
ec{iCase} = repmat(tc, nNrn, 1)
elww{iCase} = maxCoupling * ones(nNrn, 1);

iCase = 2;
ec{iCase} = allColors(round(linspace(1,120, 6)), :);
elww{iCase} = maxCoupling * ones(nNrn, 1);

iCase = 3;
ec{iCase}(1,:) = allColors(1, :);
ec{iCase}(4,:) = allColors(1, :);

ec{iCase}(3,:) = allColors(120, :);
ec{iCase}(2,:) = allColors(120, :);

ec{iCase}(5,:) = allColors(240, :);
ec{iCase}(6,:) = allColors(240, :);

elww{iCase} = maxCoupling * ones(nNrn, 1);


iCase = 4;
tc = allColors(180, :);
ec{iCase} = repmat(tc, nNrn, 1)
elww{iCase} = maxCoupling * zeros(nNrn, 1);


for iCase = 1 : nCase
    subplot2(nR,nC, iCase+1, 1);    
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
    illustrateToyModel1(ec{iCase}, elww{iCase})
    if iCase == nCase
        cbh = colorbar('southoutside');
        title(cbh, 'Phase lag [deg]')
        set(cbh, 'fontsize', 9)
        colormap(gca, 'hsv');
        displaceFigureStuff(cbh, [.017 -.08 NaN NaN]);
        caxis([0 360])
        set(cbh, 'xtick', [0 180 360])
        % set(cbh, 'Position', [0.1632 0.1003 0.1215 0.0065])
    end
end

%% subplot labels
vc.f1.sblfz = 14;

vcrd.r1y = 0.9341;
vcrd.r2y = 0.7751;
vcrd.r3y = 0.5841;
vcrd.r4y = 0.4141;
vcrd.r5y = 0.2341;

vcrd.c1x = .08;
vcrd.c2x = 0.71;

vcrd.all = [...
    vcrd.c1x vcrd.r1y; ... % A
    vcrd.c2x vcrd.r1y; ... % B
    vcrd.c1x vcrd.r2y; ... % C
    vcrd.c1x vcrd.r3y; ... % D
    vcrd.c1x vcrd.r4y; ... % E
    vcrd.c1x vcrd.r5y      % F
           ];

vcrd.labels = {'A', 'B', 'C', 'D', 'E', 'F'};


for k = 1 : size(vcrd.all, 1)
    annotation('textbox',[vcrd.all(k,1) vcrd.all(k,2)  1 0.035], ...
               'String', vcrd.labels{k}, ...
               'FontSize', vc.f1.sblfz, ... 
               'LineWidth',1, ...
               'EdgeColor', 1*[1 1 1], ...
               'BackgroundColor', [1 1 1], ...
               'FitBoxToText','on')
end

