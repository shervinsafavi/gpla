run v3_commonStuff.m

%%
signalParams.dt         = 1/signalParams.SF;               % length of time bin (in s)
signalParams.t          = 0 : signalParams.dt : signalParams.signalLength - signalParams.dt;
signalParams.nSample    = numel(signalParams.t);


% allOscFreq = [11  12 13 14 15];

 %%
% % this is useless
% noiseLFP = ...
%     globalDynamicsParams.whiteNoise_sigma ...
%     * randn(signalParams.nCh, signalParams.nSample, signalParams.nTr);

nOscComp = numel(globalDynamicsParams.oscFreq);
% amps = globalDynamicsParams.oscComps;
% oscWeights = amps / sum(amps);
% oscFreqs = globalDynamicsParams.oscFreq;

for iOscComp = 1 : nOscComp
    cmplx_oscComps(iOscComp, :) =  ...
        exp(1i * (2*pi* globalDynamicsParams.oscFreq(iOscComp) * signalParams.t));
end

%%
% clear signalParams