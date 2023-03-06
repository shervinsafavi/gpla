function [lfpLikeSig, spikeTrains, varargout]= ...
    smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, couplingParams, signalParams, varargin)
% [lfpLikeSig, spikeTrains]= smlt_transLockedSpkLfpPairs(globalDynamicsParams, spikeTrainParams, signalParams)
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
% 1b        .nCycl: scalar
%            TEXT
% 1c        .syncSigProportion: scalar
%            TEXT
% 1d        .lfpPhaseNoise_kappa: scalar
%            TEXT
% (1?)        .lfpPhaseShift: scalar, vector or a matrix
%            TEXT
% 1e        .whiteNoise_sigma: scalar
%            TEXT
% 1f        .lockingStrength_kappa: scalar, vector or a matrix
%            parameter determining the concentration 
%            (similar to what we have in a von Mises distribution)
% (1g)      .lockingPhase: scalar, vector or a matrix
%            phase (in radian) where spikes are locked in    
% 1     spikeTrainParams: structure
%       contain features of the spike train(s) and should include the following fields:
% 1a        .avefiringRate: scalar, vector or a matrix
%            average spiking rate over time in the transient epochs                  
% 2     signalParams: structure
%       contain signal params and should include the following fields:
% 2a        .signalLength: scalar
%            length of your LFP time series and spike train in second
% 2b        .SF: scalar
%            sampling frequency
% 2c        .nTrial: scalar
%            number of trials
% 2d        .nCh: scalar
%            number of channels
% 2e        .nUnit: scalar
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
%   creation: 2018-06-03 by SS (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ YYYY-MM-DD TEXT 
% ------
% see also 

%% handle optional variables
nOpVar = 0; % number of optional variable
nOpVar = nOpVar + 1; opVar.eventTrain   = []; defaultValues{nOpVar} = NaN;
 
opVar = handleVarargin(varargin, opVar, defaultValues);


%% some basics 
signalParams.dt         = 1/signalParams.SF;               % length of time bin (in s)
signalParams.t          = 0 : signalParams.dt : signalParams.signalLength - signalParams.dt;
signalParams.nSample    = numel(signalParams.t);

%% generate continues (un-)coupled-spike
PLspikeTrainParams = spikeTrainParams; % PL = phase-locked
PLspikeTrainParams.kappa        = couplingParams.lockingStrength_kappa;
PLspikeTrainParams.lockingFreq  = globalDynamicsParams.oscFreq;
PLspikeTrainParams.lockingPhase = couplingParams.lockingPhase;

PLspikeTrains  = gnrt_phaseLockedSpikeTrains(PLspikeTrainParams, signalParams); 

% % some triying to see if mult-locking will effect stuff
% warning('funny stuff here')
% tmp1  = gnrt_phaseLockedSpikeTrains(PLspikeTrainParams, signalParams); 
% tmpPLspikeTrainParams = PLspikeTrainParams;
% tmpPLspikeTrainParams.lockingFreq = PLspikeTrainParams.lockingFreq - 10.04;
% tmpPLspikeTrainParams.kappa = PLspikeTrainParams.kappa * 10;
% tmp2  = gnrt_phaseLockedSpikeTrains(tmpPLspikeTrainParams, signalParams); 
% % tmpPLspikeTrainParams.lockingFreq = PLspikeTrainParams.lockingFreq - .05;
% % tmp3  = gnrt_phaseLockedSpikeTrains(tmpPLspikeTrainParams, signalParams); 
% PLspikeTrains  = tmp1 + tmp2;
% PLspikeTrains  = tmp1 + tmp2 + tmp3;

%% generate continues osc LFP
% cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));
% oscLFP: repmat cmplx_oscComp
% tmpMat = repmat(cmplx_oscComp, signalParams.nCh, signalParams.nTr);
% noiseLFP: randn
noiseLFP = ...
    globalDynamicsParams.whiteNoise_sigma ...
    * randn(signalParams.nCh, signalParams.nSample, signalParams.nTr);

% cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));

if isfield(globalDynamicsParams, 'oscComps')  
    warning('there is more comps')
    % there are mult single freq comp
    oscComps = globalDynamicsParams.oscComps;
    nOscComp = size(oscComps, 1);
    
    % scale the amps
    tmpamps = oscComps(:,2); % second col cotains the ams, first freqs
    tmposcFreqs = oscComps(:,1);
    amps  = [1; tmpamps];
    oscWeights = amps / sum(amps);
    
    %    
    tmpfreq = oscComps(:,1); % second col cotains the ams, first freqs
    oscFreqs = [globalDynamicsParams.oscFreq; tmpfreq];
    
    % non-main freqs
    for iOscComp = 1 : nOscComp + 1
        cmplx_oscComps(iOscComp, :) = oscWeights(iOscComp) ...
            * exp(1i * (2*pi* oscFreqs(iOscComp) * signalParams.t));
    end
% need to be investigated     
    cmplx_oscComp = sum(cmplx_oscComps, 1);
% not sure if this phase noise make any sense anymore    
    cmplx_noisyOsc_raw = repmat(cmplx_oscComp, signalParams.nCh, 1, signalParams.nTr) ...
        .* exp(1i * randvm(globalDynamicsParams.lfpPhaseNoise_kappa, size(noiseLFP))) ;   
else
    % there os a single freq comp
    cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));
    cmplx_noisyOsc_raw = repmat(cmplx_oscComp, signalParams.nCh, 1, signalParams.nTr) ...
        .* exp(1i * randvm(globalDynamicsParams.lfpPhaseNoise_kappa, size(noiseLFP)));
end



if ~isfield(globalDynamicsParams, 'lfpPhaseShift')  
%     globalDynamicsParams.lfpPhShift = 0;
    cmplx_noisyOsc = cmplx_noisyOsc_raw;
else
    % shift the lfp also
    cmplx_noisyOsc = cmplx_noisyOsc_raw ...
        .* exp(1i * repmat(globalDynamicsParams.lfpPhaseShift, 1, signalParams.nSample, signalParams.nTr));
end

%% generate Poisson events
% events are indicating the center of transient osc
% compute the rate needed

transientUnit_duration = (globalDynamicsParams.nCycl * (1/globalDynamicsParams.oscFreq));

if (~iscell(opVar.eventTrain) && isnan(opVar.eventTrain))    
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
else
    eventTrain = opVar.eventTrain;
end

%% generate all-zero spike
for iTr = 1 : signalParams.nTr
    spikeTrains{iTr} = sparse(zeros(signalParams.nUnit, signalParams.nSample));
end

%% replacing event windows 
transientUnit_t = 0 : signalParams.dt : transientUnit_duration-signalParams.dt;
transientDuration_nSample = numel(transientUnit_t);
varargout{1} = transientDuration_nSample;

lfpLikeSig = noiseLFP;
% for event windows replace osc LFPs with noise (maybe with a quasssian for filtering issues)
for iTr = 1 : signalParams.nTr
    nEvent = numel(eventTrain{iTr});
    for iEvent = 1 : nEvent
        tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
        tmpRange = tmpInd : tmpInd+transientDuration_nSample-1;
%         oscLFP
        gaussWinnedLFP =  real(cmplx_noisyOsc(:,tmpRange, iTr)) ...
            .* repmat(gausswin(numel(transientUnit_t))', signalParams.nCh, 1);
        lfpLikeSig(:, tmpRange, iTr) = gaussWinnedLFP;
    end
end

% for event windows replace w/ (un-)modulated FRM
for iTr = 1 : signalParams.nTr
    nEvent = numel(eventTrain{iTr});
    for iEvent = 1 : nEvent        
        tmpInd = round(eventTrain{iTr}(iEvent) * signalParams.SF);
        tmpRange = tmpInd : tmpInd+transientDuration_nSample-1;
        spikeTrains{iTr}(:, tmpRange) =  PLspikeTrains(:, tmpRange, iTr);        
    end
end


%% oldies
% cmplx_oscComp =  exp(1i * (2*pi* globalDynamicsParams.oscFreq * signalParams.t));
% 
% 
% betaVals = abs(cmplx_betaVals);
% 
% % for iU = 1 : signalParams.nSpk 
% %     tmp_cmplx_oscComp = cmplx_oscComp * exp(1i * angle(cmplx_betaVals(iU)));
% %     tmp_oscComp = (1 + real(tmp_cmplx_oscComp)) / 2;
% %     FRmodulators(iU, :) = betaVals(iU) * tmp_oscComp + (1 - betaVals(iU)) * .5;
% % %     cmplx_FRmodulators(iU, :) = abs(betaVals(iU)) * oscComp ...
% % %         + (1 - betaVals(iU)) * .5;
% % end
% 
% % tmpNoise = noise_sigma * randn(signalParams.nLfp, signalParams.nSample);
% 
% for iCh = 1 : signalParams.nLfp        
% % 	cmplx_noisyOsc = exp(1i * (2*pi* globalWaveParams.oscFreq * signalParams.t + randvm(kappa, [1 numel(oscComp)], pi/10)));
% %     noisyOsc = alphaVals(iCh) * oscComp + tmpNoise(iCh, :); 
%     cmplx_noisyOsc = abs(cmplx_alphaVals(iCh)) ...
%         * cmplx_oscComp .* exp(1i * randvm(kappa, size(cmplx_oscComp), angle(cmplx_alphaVals(iCh))));  
% 
% % 
% %     %  *** we substract mean for proper estimation of the phase with
% %     %  Hilbert
% %     tmpAnalSig = hilbert(noisyOsc - mean(noisyOsc));
% %     lfpLikeSig_phase(iCh, :) = angle(tmpAnalSig);
% % %     lfpLikeSig(iCh, :) = oscComp;    
%     tmpLfpLikeSig(iCh, :) = real(cmplx_noisyOsc);    
% end



% %% generate firing rate modulators 
% 
% transientUnit_t = 0 : signalParams.dt : transientUnit_duration-signalParams.dt;
% transientDuration_nSample = numel(transientUnit_t);
% 
% % generate modulator unit together with the gaussian window for FR modulator
% gaussWinAlpha = (transientUnit_duration/3)^-1;
% transientUnit_t = 0 : signalParams.dt : transientUnit_duration-signalParams.dt;
% % FRmodulatorUnit = ...
% %     ((1 + cos(2*pi* globalWaveParams.oscFreq * transientUnit_t)) / 2) ...
% %     .* (gausswin(numel(transientUnit_t),gaussWinAlpha))';
% 
% % by replacing zeros with a gussain windowed osc
% FRmodulators = zeros(signalParams.nSpk, signalParams.nSample);
% 
% % to contruct the FR modu units
% % oscComp = (1 + cos(2*pi* globalWaveParams.oscFreq * signalParams.t)) / 2;
% for iU = 1 : signalParams.nSpk
%     tmp_cmplx_oscComp = cmplx_oscComp * exp(1i * angle(cmplx_betaVals(iU)));
%     tmp_oscComp = (1 + real(tmp_cmplx_oscComp)) / 2;        
%     FRmodulatorsForUnits(iU, :) = betaVals(iU) * tmp_oscComp + (1 - betaVals(iU)) * .5;
% end
% 
% for iU = 1 : signalParams.nSpk
%     for iEvent = 1 : numel(eventTrain)
%         tmpInd = round(eventTrain(iEvent) * signalParams.SF);
%         tmpFRmodulatorUnit = ...
%             FRmodulatorsForUnits(iU, tmpInd:tmpInd+transientDuration_nSample-1) ...
%             .* (gausswin(numel(transientUnit_t),gaussWinAlpha))';               
%         FRmodulators(iU, tmpInd:tmpInd+transientDuration_nSample-1) = tmpFRmodulatorUnit;       
%     end
% end
% 
% lfpLikeSig = noise_sigma * randn(signalParams.nLfp, signalParams.nSample);
% for iCh = 1 : signalParams.nLfp
%     for iEvent = 1 : numel(eventTrain)
%         tmpInd = round(eventTrain(iEvent) * signalParams.SF);
%         tmpLfpUnit = tmpLfpLikeSig(iCh, tmpInd:tmpInd+transientDuration_nSample-1);
% %         tmpFRmodulatorUnit = ...
% %             FRmodulatorsForUnits(iCh, tmpInd:tmpInd+transientDuration_nSample-1) ...
% %             .* (gausswin(numel(transientUnit_t),gaussWinAlpha))';               
%         lfpLikeSig(iCh, tmpInd:tmpInd+transientDuration_nSample-1) = tmpLfpUnit;       
%     end
%     tmpSig =  lfpLikeSig(iCh,:);
%     tmpAnalSig = hilbert(tmpSig - mean(tmpSig));
%     lfpLikeSig_phase(iCh, :) = angle(tmpAnalSig);
%     
% end
% 
% 
% %% generate modulated spikes
% for iU = 1 : signalParams.nSpk
%     spikeTrains = gnrt_inhomogeneousPoissonSpikeTrain(scaleFactor*FRmodulators, signalParams.T, signalParams.SF);
% end

%%
% clf
% tmpRange = 1300:2000; 
% plot_eventRaster(spikeTrains(1:20, tmpRange), {'color','k'})
% hold on; 
% plot_eventRaster(spikeTrains(1:5, tmpRange), {'color','r'})
% grid on