%% 
% this fig demonstrate some analysis related to statistics of GPLA
% singular values 

%%
clear all
ignit_gpla

run('../vizConventions')

%%
addpath(fullfile(dirPath_gpla.src, 'exptools'))

%%
fs = 10;
lw = 2;

%%
nRel = 2;

allOscFreq = [11  12 13 14 15];

%%
globalDynamicsParams = struct ...
    (...
        'oscFreq',  allOscFreq, ...  % Hz   % oscillation frequency for LFP and phase-locking (will be assigned later)
        'lfpPhaseNoise_kappa', ...   % concentration parameter $\kappa$ for noise on LFP phase
        10, ...
        'whiteNoise_sigma', ...      % variance for white noise for non-oscillatory periods
        0 ...  
        );
globalDynamicsParams.nFreqComp = numel(globalDynamicsParams.oscFreq);

%
signalParams = struct ...
    ( ...
        'nCh',      NaN, ...          % number of channel with LFP
        'nUnit',    NaN, ...          % number of spiking units
        'SF',       1e3, ... % Hz    % sampling frequency
        'nTr',      10, ...          % number of trials
        'signalLength', ...  % sec   % duration of the signals
        110 ...
        );


%
noc = numel(globalDynamicsParams.oscFreq);
% *** should be removed
% oscWeights are useless, given we have mixing matrix here, for
% sake syntax compatibilty  we keep them
globalDynamicsParams.oscComps = ones(1, noc);


%%
% this is the numbers we used for LFP not spiking unit
% unitNums = 2 .^ (2:11);
unitNums = (1:10)*100;
nUnitNum = numel(unitNums);

clear tv
tv = 10 .^ -(0 : 4)
couplingStrengths = sort([0 tv tv(1:end-1)/2])
% couplingStrengths = 10 .^ (-1*(5:-1:0) ) ; %
nCouplingStrength = numel(couplingStrengths);


%% filter
filterOrder = 2;
nIteCent = 0;
halfFilterWidth = 5;
iFreq = 1; % which freq of lfp osc is being used for GPLA and
           % spike-lfp analysis
           % Adjust parameters related to frequencies:
freqCenter = globalDynamicsParams.oscFreq(iFreq);

%% stats
statTestInfo.testType = 'RMT-based';
statTestInfo.SVspectrumStatsType = 'RMT-heuristic';

sameElecCheckInfo = [];

%%
%% detection as a function of coupling and dimention
caseName = 'correctPosInves'; 

% cer = []; % clusterer: define the population structure
ncc = 1;  % we only play with one case of cell clusterting
% nclu = globalDynamicsParams.nFreqComp; % number of clusters (because we deal with 6 freq)
                                       % cer = eye(nclu);
                                       % cer((3:6), (3:6)) = 0;

nclu = globalDynamicsParams.nFreqComp;

% firing rate is specified for all population
aveFR = 20;

%
iCase = 1; % is useless here

% ~ allocate 
summaryStat.(caseName) = NaN(nUnitNum, nCouplingStrength, nRel);

for iun = 1 
    signalParams.nCh = unitNums(iun);
    signalParams.nUnit = unitNums(iun);

    % coupling strength is specified later for the population
    couplingParams(iCase) = struct ...
        (...
            'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
            NaN, ...
            'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
            2*pi * rand(signalParams.nUnit, 1) ... % rad  
            );

    % for icc = 1 : ncc % loop is not necessary here
    % globalDynamicsParams.oscFreq = allOscFreq;
    
    % mixing matrix
    baseLineCoef = 0.1; % 0 correspond to no mixing
    globalDynamicsParams.mixingMatrix = kron(eye(noc), (1 - baseLineCoef)*ones(signalParams.nCh/noc,1)) + baseLineCoef;

    baseLineCoef = 0;
    spikeTrainParams.avefiringRate = kron(eye(noc), aveFR * (1 - baseLineCoef)*ones(signalParams.nUnit/noc,1)) + baseLineCoef;

    for ic = 1 
        couplingStrength = couplingStrengths(ic);
        % coupling strength is specified for all population
        couplingParams(iCase).lockingStrength_kappa = zeros(signalParams.nUnit, globalDynamicsParams.nFreqComp);
        couplingParams(iCase).lockingStrength_kappa(:, 1) = couplingStrength * ones(signalParams.nUnit, 1);


        for iRel = 1
            iRel
            %% simulation 
            [lfpLikeSig.(caseName), spikeTrains.(caseName), spikeTrains_xPh.(caseName)] = ...
                tsmlt_sustMixedLockedSpkLfpPairs_multFreqPureCoupling5(...
                    globalDynamicsParams, spikeTrainParams, couplingParams, ...
                    signalParams);
            
            %% filter LFP 
            [~, analLfp.(caseName)] = tpp_filt_recenter(lfpLikeSig.(caseName), freqCenter, halfFilterWidth, ...
                                                        signalParams.SF, filterOrder, nIteCent);

            %% GPLA
            icc = 1; % usless index, should be removed
            caseNames{1} = caseName;

            [~,~,~,~, allStats.(caseName){iun, iRel}, svdOut.(caseName){iun, iRel}] = ...
                tngpla(spikeTrains_xPh.(caseName), analLfp.(caseName), [], [], [], [] , [], ...
                       statTestInfo, 'all', sameElecCheckInfo, ...
                       'nSpk-square-root', 2, 0);
        end
        
        % summaryStat.(caseName)(iun, ic, iRel) = allStats.falsePosInves{iun,iRel}.gPLV_stats.nullHypoReject;
        
    end
end                                     %

