vc = get_vizConventions();
storagePath = fullfile('.','dat');

%%
nRel = 100;

allOscFreq = [11  12 13 14 15];

mixingBaseLineCoef = .1;

%%
lw = 2.5;

%%
signalParams = struct ...
    ( ...
        'nCh',      NaN, ...          % number of channel with LFP
        'nUnit',    NaN, ...          % number of spiking units
        'SF',       1e3, ... % Hz    % sampling frequency
        'nTr',      10, ...          % number of trials
        'signalLength', ...  % sec   % duration of the signals
        11 ...
        );

% firing rate is specified for all population
aveFR = 20;

%% we need to change the parameters to have phase noise
clear globalDynamicsParams

globalDynamicsParams = struct ...
    (...
        'oscFreq',  allOscFreq, ...  % Hz   % oscillation frequency for LFP and phase-locking (will be assigned later)
        'lfpPhaseNoise_kappa', ...   % concentration parameter $\kappa$ for noise on LFP phase
        10, ...
        'whiteNoise_sigma', ...      % variance for white noise for non-oscillatory periods
        0 ...  
        );
noc = numel(globalDynamicsParams.oscFreq);
globalDynamicsParams.nFreqComp = numel(globalDynamicsParams.oscFreq);
globalDynamicsParams.oscComps = ones(1, noc);


%%
% *** should be removed
% oscWeights are useless, given we have mixing matrix here, for
% sake syntax compatibilty  we keep them
% globalDynamicsParams.oscComps = ones(1, noc);

%%
% this is the numbers we used for both LFP and spiking units
% unitNums = 2 .^ (2:11);
unitNums = (1:10)*100;
nUnitNum = numel(unitNums);

clear tv
tv = 10 .^ -(0 : 4)
couplingStrengths = sort([0 tv tv(1:end-1)/2])
% couplingStrengths = 10 .^ (-1*(5:-1:0) ) ; %
nCouplingStrength = numel(couplingStrengths);

%%
popNums = 1 : 10;
nPopNum = numel(popNums);


%% filter
filterOrder = 2;
nIteCent = 5;
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
nR = 3;
nC = 4;
