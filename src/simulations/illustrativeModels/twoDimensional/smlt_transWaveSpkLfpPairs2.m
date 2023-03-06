function [lfpLikeSig, spikeTrains, map, varargout]= ...
    smlt_transWaveSpkLfpPairs2(globalDynamicsParams, spikeTrainParams, ...
    pdeParams, sigmfParams)


%%
SF = globalDynamicsParams.oscFreq * globalDynamicsParams.stimPer;

%% generate waves 
% generate a sequence of transients start times
transientParams = globalDynamicsParams;
transientParams.stimSequence = ...
    transientParams.waveFormType * ones(1, transientParams.nTransients);

[stimVec, ~] = generWavStim(transientParams);
varargout{2} = stimVec;

waveData = sim_WaveEqu(stimVec,'edp_params',pdeParams);

%%
% We reshape the wave data and generate it corresponding map
map = nan(size(waveData, 1), size(waveData, 2));
unitCounter = 0;
for iX = 1 : size(waveData, 1)
    for iY = 1 : size(waveData, 2)
        unitCounter = unitCounter + 1;
        % create the map
        map(iY, iX) = unitCounter;        
    end
end

nUnit   = unitCounter; clear unitCounter 
nSample = size(waveData, 3);
nTr     = 1;

%% LFP-like signals
% waveData_reshaped = reshape(waveData, nUnit, nSample);

for it = 1 : nSample
% for it = 1000% : nSample
    tm = waveData(:,:, it);
    lfpLikeSig_inGridShape(:,:, it) = imgaussfilt(tm, 2) + globalDynamicsParams.whiteNoise_sigma * randn(size(tm));
end

lfpLikeSig = reshape(lfpLikeSig_inGridShape, nUnit, nSample);

% subplot 121
% imagesc(tm)
% colorbar
% axis square
% subplot 122
% imagesc(tms)
% colorbar
% axis square


%%

% lfpLikeSig = waveData_reshaped ...
%     + globalDynamicsParams.whiteNoise_sigma * randn(size(waveData_reshaped));

%% spike trains
waveMoulator = nan(nUnit, nSample);
waveData_reshaped = reshape(waveData, nUnit, nSample);
varargout{4} = waveData_reshaped;


% parameters for the transformation with the sigmoid function:
% sigmf(X, [A, C]) = 1./(1 + EXP(-A*(X-C)))

% sigmfParams.A            = 2;
% sigmfParams.C            = 0;
% sigmfParams.sharpenerExp = 40;  % heuristic for sharpening the function

for iUnit = 1 : nUnit
    waveMoulator(iUnit, :) = ...
        (sigmf(waveData_reshaped(iUnit, :), ...
        [sigmfParams.A sigmfParams.C])) .^ sigmfParams.sharpenerExp;       
end

if numel(spikeTrainParams.betaVal) == 1
    
    betaVal = spikeTrainParams.betaVal;
    FRmoulator = spikeTrainParams.scaleParam ...
        * (betaVal * waveMoulator + (1 - betaVal) * spikeTrainParams.baselineFR);
    spikeTrains_binMat = ...
        gnrt_inhomogeneousPoissonSpikeTrains(FRmoulator, nSample/SF, SF, nTr);
    
    varargout{3} = FRmoulator; 
    
    %%
    for iTr = 1 : nTr
        % org spike trains
        spikeTrains{iTr} = sparse(spikeTrains_binMat(:,:, iTr));
    end
else    
    betaVals = spikeTrainParams.betaVal;
    nBeta = numel(betaVals);
    tmpFR{nBeta} = nan(nUnit, nSample);
    for iBeta = 1 : nBeta
        tmpFR{iBeta} = nan(nUnit, nSample);
        betaVal = betaVals(iBeta);
        FRmoulator = spikeTrainParams.scaleParam ...
            * (betaVal * waveMoulator + (1 - betaVal) * spikeTrainParams.baselineFR);
        spikeTrains_binMat = ...
            gnrt_inhomogeneousPoissonSpikeTrains(FRmoulator, nSample/SF, SF, nTr);
        tmpFR{iBeta} = FRmoulator; 
        
        for iTr = 1 : nTr
            spikeTrains{iBeta}{iTr} = sparse(spikeTrains_binMat(:,:, iTr));
        end
    end
    varargout{3} = tmpFR;
end


%%
signalParams.nUnit      = nUnit;
signalParams.nSample    = nSample;
signalParams.nTr        = nTr;
signalParams.SF         = SF;

varargout{1} = signalParams;