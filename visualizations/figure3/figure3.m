%% Figure 3: Comparing GPLA and PLA
% This figure demonstrate the advantegious of GPLA to PLA (when both are
% applicabale)

%% initiatation 
% (e.g add necessary packages and functions)clear all
clear all
clf
pds = ignit();

%% prepare the figure

% visualization conventions
vc = get_vizConventions();

fh = gcf;

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f9.w;
fh.Position(4)  = vc.f9.h;

% number of rows and columns in the subplot
nR = 4; 
nC = 4;

%% using different number of trials
% computation will take long, so the resulting mat file is load

% one can run the scrip to do the computation though
% run v5_plvEstCmpr_trial.m
load(fullfile('dat', 'v5_plvEstCmpr_trial.mat'))
tos = -.035; % temp offset to move last 2 subplot down
vc = get_vizConventions();
nR = 4; nC = 4;

%%
[~, sbh] = subplot2d(nR,nC, [3 4], [1 2])
displaceFigureStuff(sbh, [-.005 tos NaN NaN])
run plvEstCmpr_trial_v5_plot_v2.m

%% using different firing rates
% computation will take long, so the resulting mat file is load

% one can run the scrip to do the computation though
% run v5_plvEstCmpr_fr.m
load(fullfile('dat','v5_plvEstCmpr_fr'));
vc = get_vizConventions();

%% first run the example stuff and then the main plot of gpla vs
%% pla comparaiosn
%% noisy vs noiseless oscillation 
% this is for some adjusment reasons
subplot2(nR,nC, 2, 3);

T = 4;
dt = .001;
t = 0:dt:T;
t(end) = [];

iCh = 1;

ksigma = 7;
sigma = lfpAmpNoise_sigmas(ksigma);
nlfp = 1;
tlfp = (ones(nlfp, 1) * exp(1i * 2 * pi * t)) + sigma*(randn(nlfp, length(t))+1i*randn(nlfp,length(t)));
plot(real(tlfp(iCh, :)), 'k')
hold on

tlfp = (ones(nlfp, 1) * exp(1i * 2 * pi * t));
plot(real(tlfp(iCh, :)), 'b', 'linewidth', 3)
hold on

txr = numel(tlfp); % temp x range
xlim([0 txr])
set(gca, 'xtick', [0 txr/2 txr])
tyr = 25;
ylim(tyr * [-1 1])
set(gca, 'ytick', [-tyr 0  tyr])
xlabel('Time [ms]')

ylabel('Amplitude [a.u.]')

[lh, icons] = legend('Noisy', 'Original', 'location','northeast')
set(lh, 'box', 'off')
displaceFigureStuff(lh, [0.015 .04 0 0])


bof = 0.3; % backward offset
for k = 3 : 2: 6
    tmpOrigPos = icons(k).XData(2);
    icons(k).XData(2) = tmpOrigPos - bof; 
end

for k = 1 : 2
    icons(k).Position(1) = 0.3;
end

set(gca, 'fontsize', vc.f9.gfs)
set_ticksOutward1(.05)
box off

%% example coupling matrix
% this should come after simulation related to checking for GPLA vs PLA plotting
subplot2(nR,nC, 2, 4)
imagesc(abs(rawSvdStuff.couplingMatrix));
set(gca, 'xtick', (10:20:signalParams.nUnit))
set(gca, 'ytick', (10:10:signalParams.nCh))
colormap(gca, 'copper')
c = colorbar; title(c, 'PLV'); 
displaceFigureStuff(c, [.06 NaN NaN NaN])
ylabel('LFP ID')
xlabel('Spiking unit ID')

set(gca, 'fontsize', vc.f9.gfs)


%
[~, sbh] = subplot2d(nR,nC, [3 4], [3 4]);
displaceFigureStuff(sbh, [NaN tos NaN NaN])
run plvEstCmpr_fr_v5_plot_v2.m

%% compare PLV and gPLV
% run comparePLVandGPLA_v1.m
load(fullfile('dat','comparePLVandGPLA_v1'));
vc = get_vizConventions();
run comparePLVandGPLA_v0_plot_v1.m

%%
% We simulate Poisson spike trains transiently phase-locked  to an
% oscilation in a prticular freqency band.

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
freqBand = globalDynamicsParams.oscFreq + [-2.5 +2.5];
wn = freqBand / (signalParams.SF/2);
filterOrder = 2;
[b, a] = butter(filterOrder, wn, 'bandpass');

%% define stat info
statTestInfo = [];
sameElecCheckInfo = [];

%% viz conve
iTr_forViz = 1;
iEvent_forViz = 2;
tmpOffSet = 100;
transDur = transientUnit_duration * signalParams.SF;


%%
nCase = 1;

%% toy simulation for demo

% define parameters
iCase = 1;
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
        5 , ...
        'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
        phaseShifter(:) ... % rad  
        );

caseType = 'gplaDemo';
iSV = 1;

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

[lfpVec.(caseType)(:,iCase), spkVec.(caseType)(:,iCase), gPLV.(caseType)(iCase), ...
 ~, allStats.(caseType)(iCase), rawSvdStuff.(caseType)(iCase)] = ...
    tngpla(spikeTrains.(caseType){iCase}, lfpLikePhase.(caseType){iCase}, [], [], [], [] , [], ...
           statTestInfo, iSV, sameElecCheckInfo, 'nSpk', 0);


%% some scripts only needed for visualization with MATLAB 
% (e.g. arranging subplots etc)
iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpRng = tmpInd-tmpOffSet : tmpInd+transDur+tmpOffSet ;

iTr  = iTr_forViz;
iEvent = iEvent_forViz;
tmpOffSet = 70;

tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
tmpMid = round(eventTrain{iTr}(iEvent) * signalParams.SF ...
               + transDur / 2);
tmpRng = tmpMid-tmpOffSet : tmpMid+tmpOffSet ;




subplot2(nR,nC, 1, [1 2])
ttr = tmpInd-100:tmpInd+3000;
tlfp = lfpLikeSig.(caseType){iCase}(iCh,ttr , iTr);

plot(ttr, tlfp, ...
     'color','k', 'linewidth', 1); %colorbar


axis tight

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))


tyr = 1.3 * [-1 1]; % temp Y range
ylim(tyr);

yslos = .1;

ystos = .11;
stoslw = 1.5;
slosfs = 8.5;

plotScaleBar1(ttr, tlfp, ...
              '500 ms', NaN, 500, NaN, ...
              yslos, NaN, ystos,  NaN, ...
              stoslw, slosfs); 


axis off
box off



tlfp = lfpLikeSig.(caseType){iCase}(1,tmpRng, iTr);
rectangle(gca, 'Position',[tmpRng(1) -1.1 numel(tmpRng) 2.2], ...
          'EdgeColor', 'b', 'linewidth', 1.2);

yslos = -10.9;
ystos = -.905;

subplot2(nR,nC, 1, 3)


trange = [(1:6); (7:12); (13:18)];
allK = [1 2 3];

clusterCols = eye(3);
clusterCols(3,:) = [1 0 1];

for k = 1 : 3
    tspt = full(spikeTrains.(caseType){iCase}{iTr}(:,tmpRng));
    excludeK = setdiff(allK, k);
    for kexc = 1 : numel(excludeK)
        tspt(trange(excludeK(kexc), :),:) = 0;
    end
    plot_EveConData(tlfp, tspt, ...
                    {'color','k', 'linewidth',vc.lw_lfp}, {'color',clusterCols(k, :), ...
                        'linewidth',vc.lw_spk})        

end

plotScaleBar1([1 numel(tlfp)], tlfp, ...
              '50 ms', NaN, 50, NaN, ...
              yslos, NaN, ystos,  NaN, ...
              stoslw, slosfs); 


set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%%
subplot2(nR,nC, 1,4);
for k = 1 : 3
    tsv = spkVec.(caseType)(trange(k,:),iCase);
    polarplot(angle(tsv), abs(tsv), ...
              'color',clusterCols(k,:), 'LineStyle','none', ...
              'marker','.');
    hold on
end

thetaticks(0:45:315)
tn = .5;
rlim([0 tn])
rticks([0 tn/2 tn]);

ax = gca;
ax.LineWidth = 1.5;
rruler = ax.RAxis;
ax.RAxisLocation = 45/2;
set(gca, 'fontsize', vc.f9.gfs)

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%% subplot labels
vc.f11.sblfz = 14;

vcrd.r1y = 0.9541;
vcrd.r2y = 0.735;
vcrd.r3y = 0.505;

vcrd.c1x = 0.05;
vcrd.c2x = .286;
vcrd.c3x = .50;
vcrd.c4x = .705;

vcrd.all = [...
    vcrd.c1x vcrd.r1y; ... % A
    vcrd.c3x vcrd.r1y; ... % B
    vcrd.c4x vcrd.r1y; ... % C
    vcrd.c1x vcrd.r2y; ... % D
    vcrd.c2x vcrd.r2y; ... % E
    vcrd.c3x vcrd.r2y; ... % F
    vcrd.c4x vcrd.r2y; ... % G
    vcrd.c1x vcrd.r3y; ... % H
    vcrd.c3x vcrd.r3y      % I
           ];

vcrd.labels = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'};


for k = 1 : size(vcrd.all, 1)
    annotation('textbox',[vcrd.all(k,1) vcrd.all(k,2)  1 0.035], ...
               'String', vcrd.labels{k}, ...
               'FontSize', vc.f11.sblfz, ... 
               'LineWidth',1, ...
               'EdgeColor', 1*[1 1 1], ...
               'BackgroundColor', [1 1 1], ...
               'FitBoxToText','on')
end


