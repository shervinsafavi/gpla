function [spikeTrain, varargout]  = gnrt_phaseLockedSpikeTrains(spikeTrainParams, signalParams, varargin)
% [spikeTrain, PLV]  = gnrt_phaseLockedSpikeTrains(spikeTrainParams, spikeTrainDuraion, SF, nTr, nUnit)
% Generate a phase-locked Poisson spike train with the given average firing
% rate, dispersion (defined by kappa) and locking phase. The code (and notation) 
% is based on: 
% 
% [1] Ashida, G., Wagner, H. & Carr, C. E. 
% Processing of Phase-Locked Spikes and Periodic Signals. 
% in Analysis of Parallel Spike Trains 59–74 (Springer, Boston, MA, 2010). 
% doi:10.1007/978-1-4419-5675-0_4
% 
%
% EXAMPLE:
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     spikeTrainParams: structure
%       contain features of the spike train(s) and should include the following fields:
% 1a        .avefiringRate: scalar, vector or a matrix
%            average spiking rate over time                   
% 1b        .kappa: scalar, vector or a matrix
%            parameter determining the concentration 
%            (similar to what we have in a von Mises distribution)
% 1c        .lockingFreq: scalar, vector or a matrix
%            frequency (in hertz) with which spikes are oscilating
% (1d)      .lockingPhase: scalar, vector or a matrix
%            phase (in radian) where spikes are locked in    
% 2     signalParams: structure
%       contain signal params and should include the following fields:
% 2a        .signalLength: scalar
%            length of your spike train in second
% 2b        .SF: scalar
%            sampling frequency
% 2c        .nTrial: scalar
%            number of trials
% 2d        .nUnit: scalar
%            number of units(i.e. generate Poisson spike train for a population
% 
% Output:
% 1     spikeTrain: binary vector
% (2)   theoPLV: 
%
% ------
% potential improvments:
% (1) TEXT
% ------
% Code Info:
%   creation: 2018-05-18 by SS (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ YYYY-MM-DD TEXT 
% ------
% see also gnrt_homogeneousPoissonSpikeTrains
% gnrt_inhomogeneousPoissonSpikeTrains

%% Handle optional inputs (varargin):
% nOpVar = 0; % number of optional variable
% nOpVar = nOpVar + 1; optionalVariables.nTr	 = []; defaultValues{nOpVar} = size(spikeTrainParams.avefiringRate, 2); 
% nOpVar = nOpVar + 1; optionalVariables.nUnit = []; defaultValues{nOpVar} = size(spikeTrainParams.avefiringRate, 1); 
% 
% optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

if ~isfield(spikeTrainParams, 'lockingPhase')  
    spikeTrainParams.lockingPhase = 0;
end

%% 
structunpack(spikeTrainParams);

%% 
% signal params
nBin   = signalParams.signalLength * signalParams.SF;
dt     = 1 / signalParams.SF;               % length of time bin (in s)
nUnit  = signalParams.nUnit;
nTr    = signalParams.nTr;


% spike train params 
a_raw       = lockingPhase;
k_raw       = kappa;
avFR_raw    = avefiringRate;        % note by \bar{R} in cited reference [1]
f_raw       = lockingFreq;

% signalParams.t          = 0 : signalParams.dt : signalParams.T - signalParams.dt;
% signalParams.nSample    = numel(signalParams.t);

t = linspace(0, signalParams.signalLength, nBin);

%% firing rate modulator 
% function handels
p   = @(x, a, k, I_0k) exp(k * cos(x - a)) / (2 * pi * I_0k);
I_0 = @(k) integral(@(x) exp(k * cos(x)), -pi, pi) / (2*pi);
theoPLVcalculator = @(k, I_0k) integral(@(x) exp(k * cos(x)) .* cos(x), -pi, pi) / (2*pi * I_0k);

% if any of signal params varies across units
if (size(a_raw, 1) > 1 || size(k_raw, 1) > 1  || size(avFR_raw, 1) > 1  || size(f_raw, 1)  > 1)
    if (size(a_raw, 2) > 1 || size(k_raw, 2) > 1  || size(avFR_raw, 2) > 1  || size(f_raw, 2)  > 1)        
        % repeat those params which are not varying  accross and/or units
        a_all    = repeatParams(a_raw, nUnit, nTr);
        k_all    = repeatParams(k_raw, nUnit, nTr);
        avFR_all = repeatParams(avFR_raw, nUnit, nTr);
        f_all    = repeatParams(k_raw, nUnit, nTr);
        
        for iUnit = 1 : nUnit
            for iTr = 1 : nTr                
                a = a_all(iUnit, iTr); k = k_all(iUnit, iTr); avFR = avFR_all(iUnit, iTr); f = f_all(iUnit, iTr);
                
                I_0k = I_0(k);
                p_ak = p(2*pi * f * t, a, k, I_0k);
                FRmodulators(iUnit, :, iTr) = 2*pi * avFR * p_ak;
%                 FRmodulator(iUnit, :) = 2*pi * avFR * p_ak;
                theoPLV(iUnit, iTr) = theoPLVcalculator(k, I_0k);
                
%                 FRmodulators = repmat(FRmodulator, 1,1, nTr);
            end
        end
    else
        % repeat those params which are not varying  accross units
        if size(a_raw, 1)    == 1, a_all    = repmat(a_raw, nUnit, 1);    else, a_all    = a_raw;    end
        if size(k_raw, 1)    == 1, k_all    = repmat(k_raw, nUnit, 1);    else, k_all    = k_raw;    end
        if size(avFR_raw, 1) == 1, avFR_all = repmat(avFR_raw, nUnit, 1); else, avFR_all = avFR_raw; end
        if size(f_raw, 1)    == 1, f_all    = repmat(f_raw, nUnit, 1);    else, f_all    = f_raw;    end
        
        for iUnit = 1 : nUnit                       
            a = a_all(iUnit); k = k_all(iUnit); avFR = avFR_all(iUnit); f = f_all(iUnit);
            
            I_0k = I_0(k);
            p_ak = p(2*pi * f * t, a, k, I_0k);
            FRmodulator(iUnit, :) = 2*pi * avFR * p_ak;
            theoPLV(iUnit) = theoPLVcalculator(k, I_0k);
        end
        FRmodulators = repmat(FRmodulator, 1,1, nTr);
    end
elseif (size(a_raw, 1) == 1 && size(k_raw, 1) == 1  && size(avFR_raw, 1) == 1  && size(f_raw, 1)  == 1)
    % compute on FRmodulator and repeat for all
    a = a_raw; k = k_raw; avFR = avFR_raw; f = f_raw;
    I_0k = I_0(k);
    p_ak = p(2*pi * f * t, a, k, I_0k);
    FRmodulator = 2*pi * avFR * p_ak;
    FRmodulators = repmat(FRmodulator, nUnit,1, nTr);
    theoPLV = theoPLVcalculator(k, I_0k);
end

% I_0k = I_0(k);
% p_ak = p(2*pi * f * t, a, I_0k);
% FRmodulator = 2*pi * avFR * p_ak;
spikeTrain = gnrt_inhomogeneousPoissonSpikeTrains(FRmodulators, signalParams.signalLength, signalParams.SF);

% theoPLV = theoPLVcalculator(k, I_0k);
varargout{1} = theoPLV;
varargout{2} = FRmodulators;


%%
% %% handle optional inputs (varargin):
% nOpVar = 0; % number of optional variable
% nOpVar = nOpVar + 1; optionalVariables.lockingPhase      = []; defaultValues{nOpVar} = 0;
% 
% optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);
% 
% %% 
% SF      = signalParams.SF;
% nBin    = signalParams.duration * SF;
% a       = optionalVariables.lockingPhase;
% k       = kappa;
% avFR    = avefiringRate;        % note by \bar{R} in cited reference [1]
% f       = lockingFreq;
% dt      = 1 / SF;               % length of time bin (in s)
% 
% % signalParams.t          = 0 : signalParams.dt : signalParams.T - signalParams.dt;
% % signalParams.nSample    = numel(signalParams.t);
% 
% t = linspace(0, signalParams.duration, nBin);
% 
% %% firing rate modulator 
% p = @(x, a, I_0k) exp(k * cos(x - a)) / (2 * pi * I_0k);
% % I_0_integrand = @(x) exp(k * cos(x));
% I_0 = @(k) integral(@(x) exp(k * cos(x)), -pi, pi) / (2*pi);
% 
% I_0k = I_0(k);
% 
% 
% FRmodulator = 2*pi * avFR * p(2*pi * f * t, a, I_0k);

end

function paramArray_wRepeat = repeatParams(rawParamArray, nUnit, nTr)
param_arraySize = size(rawParamArray);
if param_arraySize == [1 1]    
        paramArray_wRepeat = repmat(rawParamArray, nUnit, nTr);
elseif param_arraySize == [nUnit 1]
        paramArray_wRepeat = repmat(rawParamArray, 1, nTr);
elseif param_arraySize == [nUnit nTr]
        paramArray_wRepeat = rawParamArray;
elseif param_arraySize == [1 nTr]
        paramArray_wRepeat = repmat(rawParamArray, nUnit, 1);
else
    error('wired')
end
% switch param_arraySize
%     case [1 1]
%         paramArray_wRepeat = repmat(rawParamArray, nUnit, nTr);
%     case [nUnit 1]
%         paramArray_wRepeat = repmat(rawParamArray, 1, nTr);
%     case [nUnit nTr]
%         paramArray_wRepeat = rawParamArray;
%     case [1 nTr]
%         paramArray_wRepeat = repmat(rawParamArray, nUnit, 1);
%     otherwise
% end
end
