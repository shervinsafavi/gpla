function [lfpLikeSig, spikeTrains, varargout]= ...
    smlt_sustLockedSpkLfpPairs_multFreqCoupling2(globalDynamicsParams, spikeTrainParams, couplingParams, signalParams)
% similar to smlt_sustLockedSpkLfpPairs_multFreqCoupling, but with
% some changes in the noise and we also have LFP MIXING
% [lfpLikeSig, spikeTrains]= smlt_sustLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, signalParams)
% TEXT TEXT TEXT
% 
%
% EXAMPLE:
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     globalDynamicsParams: structure
%       characterize the global dynamics and should include the following fields:
% 1a        .oscFreq: scalar
%            oscilation frequency (in hertz) of locked spikes and LFPs
% 1b        .lfpPhaseNoise_kappa: scalar
%            TEXT
% 1c        .lockingStrength_kappa: scalar, vector or a matrix
%            parameter determining the concentration 
%            (similar to what we have in a von Mises distribution)
% (1d)      .lockingPhase: scalar, vector or a matrix
%            phase (in radian) where spikes are locked in    
% 1     spikeTrainParams: structure
%       contain features of the spike train(s) and should include the following fields:
% 2a        .avefiringRate: scalar, vector or a matrix
%            average spiking rate over time in the transient epochs                  
% 3     signalParams: structure
%       contain signal params and should include the following fields:
% 3a        .signalLength: scalar
%            length of your LFP time series and spike train in second
% 3b        .SF: scalar
%            sampling frequency
% 3c        .nTr: scalar
%            number of trials
% 3d        .nCh: scalar
%            number of channels
% 3e        .nUnit: scalar
%            number of spiking units 
% 
% Output:
% 1     lfpLikeSig: 
% 2     spikeTrains: binnary matrix 
%
% ------
% potential improvments:
% (1) TEXT
% ------
% Code Info:
%   creation: 2018-07-10 by SS (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ YYYY-MM-DD TEXT 
% ------
% see also 

%% handle optional variables
% nOpVar = 0; % number of optional variable
% nOpVar = nOpVar + 1; opVar.eventTrain   = []; defaultValues{nOpVar} = NaN;
%  
% opVar = handleVarargin(varargin, opVar, defaultValues);



%% some basics 
signalParams.dt         = 1/signalParams.SF;               % length of time bin (in s)
signalParams.t          = 0 : signalParams.dt : signalParams.signalLength - signalParams.dt;
signalParams.nSample    = numel(signalParams.t);

%% generate continues (un-)coupled-spike
% ~ allocate memory for

% PLspikeTrains = zeros(signalParams.nUnits, signalParams.nSample, signalParams.nTr);
for ifreq = 1 : numel(globalDynamicsParams.oscFreq)
    % PLspikeTrainParams = spikeTrainParams; % PL = phase-locked
    PLspikeTrainParams.avefiringRate = spikeTrainParams.avefiringRate(:, ifreq); % PL = phase-locked
    PLspikeTrainParams.kappa        = couplingParams.lockingStrength_kappa(:,ifreq);
    PLspikeTrainParams.lockingFreq  = globalDynamicsParams.oscFreq(ifreq);
    PLspikeTrainParams.lockingPhase = couplingParams.lockingPhase;    
    % *** should be double check if thats going to be general enough
    [PLspikeTrains_allFreqSep(:,:,:, ifreq)]  = gnrt_phaseLockedSpikeTrains(PLspikeTrainParams, signalParams); %
end

PLspikeTrains  = logical(sum(PLspikeTrains_allFreqSep, 4));
spikeTrains_binary = PLspikeTrains;

spikeTrains{signalParams.nTr} = [];
for iTr = 1 : signalParams.nTr
    spikeTrains{iTr} = sparse(spikeTrains_binary(:,:, iTr));
end

%% generate continues osc LFP
% cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));
% oscLFP: repmat cmplx_oscComp
% tmpMat = repmat(cmplx_oscComp, signalParams.nCh, signalParams.nTr);
% noiseLFP: randn

% probably it wont be needed (except the size of which is used later)
noiseLFP = ...
    globalDynamicsParams.whiteNoise_sigma ...
    * randn(signalParams.nCh, signalParams.nSample, signalParams.nTr);

% cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));

% warning('there is more comps')
nOscComp = numel(globalDynamicsParams.oscComps);
amps = globalDynamicsParams.oscComps;
oscWeights = amps / sum(amps);
oscFreqs = globalDynamicsParams.oscFreq;

for iOscComp = 1 : nOscComp
    cmplx_oscComps(iOscComp, :) = oscWeights(iOscComp) ...
        * exp(1i * (2*pi* oscFreqs(iOscComp) * signalParams.t));
end
% % need to be investigated
% cmplx_oscComp = sum(cmplx_oscComps, 1);

% oscWeights are useless, given we have mixing matrix here
for iOscComp = 1 : nOscComp
    cmplx_oscComps(iOscComp, :) = oscWeights(iOscComp) ...
        * exp(1i * (2*pi* oscFreqs(iOscComp) * signalParams.t));
end

cmplx_osc_raw = globalDynamicsParams.mixingMatrix * cmplx_oscComps;


% % not sure if this phase noise make any sense anymore
% cmplx_noisyOsc_raw = repmat(cmplx_oscComp, signalParams.nCh, 1, signalParams.nTr) ...
%     .* exp(1i * randvm(globalDynamicsParams.lfpPhaseNoise_kappa, size(noiseLFP))) ;

% we add the noise as normal noise to both real and imaginary part
cmplx_noisyOsc_raw = repmat(cmplx_osc_raw, 1, 1, signalParams.nTr) ...
    + (randn(size(noiseLFP)) + 1i * randn( size(noiseLFP))) * globalDynamicsParams.whiteNoise_sigma ;


if ~isfield(globalDynamicsParams, 'lfpPhaseShift')  
    %     globalDynamicsParams.lfpPhShift = 0;
    cmplx_noisyOsc = cmplx_noisyOsc_raw;
else
    % shift the lfp also
    cmplx_noisyOsc = cmplx_noisyOsc_raw ...
        .* exp(1i * repmat(globalDynamicsParams.lfpPhaseShift, 1, signalParams.nSample, signalParams.nTr));
end



% %% generate continues osc LFP
% cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));
% % oscLFP: repmat cmplx_oscComp
% % tmpMat = repmat(cmplx_oscComp, signalParams.nCh, signalParams.nTr);
% % % noiseLFP: randn
% % noiseLFP = ...
% %     globalDynamicsParams.whiteNoise_sigma ...
% %     * randn(signalParams.nCh, signalParams.nSample, signalParams.nTr);

% cmplx_noisyOsc = repmat(cmplx_oscComp, signalParams.nCh, 1, signalParams.nTr) ...
%     .* exp(1i * randvm(globalDynamicsParams.lfpPhaseNoise_kappa, ...
%     [signalParams.nCh signalParams.nSample signalParams.nTr]));  

lfpLikeSig = real(cmplx_noisyOsc);