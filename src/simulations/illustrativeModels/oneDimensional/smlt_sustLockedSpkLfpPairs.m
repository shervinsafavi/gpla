function [lfpLikeSig, spikeTrains, varargout]= ...
    smlt_sustLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, couplingParams, signalParams)
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
PLspikeTrainParams = spikeTrainParams; % PL = phase-locked
PLspikeTrainParams.kappa        = couplingParams.lockingStrength_kappa;
PLspikeTrainParams.lockingFreq  = globalDynamicsParams.oscFreq;
PLspikeTrainParams.lockingPhase = couplingParams.lockingPhase;

[spikeTrains_binary, theoPLV] = gnrt_phaseLockedSpikeTrains(PLspikeTrainParams, signalParams); 
varargout{1} = theoPLV;

spikeTrains{signalParams.nTr} = [];
for iTr = 1 : signalParams.nTr
    spikeTrains{iTr} = sparse(spikeTrains_binary(:,:, iTr));
end

%% generate continues osc LFP
cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));
% oscLFP: repmat cmplx_oscComp
% tmpMat = repmat(cmplx_oscComp, signalParams.nCh, signalParams.nTr);
% % noiseLFP: randn
% noiseLFP = ...
%     globalDynamicsParams.whiteNoise_sigma ...
%     * randn(signalParams.nCh, signalParams.nSample, signalParams.nTr);

% *** instead of assigning NaN is better we fill the structre with
% some default vals at very begining and if is NaN forg3et about noise
if ((isfield(globalDynamicsParams, 'lfpPhaseNoise_kappa') && ...
    isnan(globalDynamicsParams.lfpPhaseNoise_kappa)))
    globalDynamicsParams = rmfield(globalDynamicsParams, ...
                                   'lfpPhaseNoise_kappa');
end



if isfield(globalDynamicsParams, 'lfpPhaseNoise_kappa')
    cmplx_noisyOsc = repmat(cmplx_oscComp, signalParams.nCh, 1, signalParams.nTr) ...
        .* exp(1i * randvm(globalDynamicsParams.lfpPhaseNoise_kappa, ...
                           [signalParams.nCh signalParams.nSample signalParams.nTr]));  
else
    cmplx_noisyOsc = repmat(cmplx_oscComp, signalParams.nCh, 1, signalParams.nTr);  
end

lfpLikeSig = real(cmplx_noisyOsc);