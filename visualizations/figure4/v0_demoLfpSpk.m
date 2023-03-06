% subplot2(nR,nC, 1, [2 3]);
% displaceFigureStuff(sbh, [r1_hos NaN NaN NaN])

%%
% clf %
tr = 1 : 700;
verticalOffset = 2.5;
slw = 1.5;
llw = 2.3;
sl = .4;

iu = 1; 
iTr = 1;
continuousData = lfpLikeSig(iu, tr, iTr);
eventData = spikeTrains{iTr}(iu, tr);

[nDim_cd, ~] = size(continuousData);
[nDim_ed, ~] = size(eventData);


%
% make the continuousData possitive if needed
if sum(continuousData < 0) > 0
    tmpShiftVal = -1 * min(continuousData); 
else
    tmpShiftVal = 0;
end
continuousData = continuousData + tmpShiftVal;
max_conData = max(continuousData);

continuousData_scaled = continuousData / max_conData;

%
plot(continuousData_scaled, 'LineWidth',llw, 'Color','k')

% hold all
tmpEventTime = find(eventData);
for iEvent = 1 : length(tmpEventTime)
    line([tmpEventTime(iEvent) tmpEventTime(iEvent)], 1.2+[0 sl], ...
         'color','r', 'LineWidth',slw);
end
hold all
%
iu = 100; 
iTr = 1;
continuousData = lfpLikeSig(iu, tr, iTr);
eventData = spikeTrains{iTr}(iu, tr);

[nDim_cd, ~] = size(continuousData);
[nDim_ed, ~] = size(eventData);


%
% make the continuousData possitive if needed
if sum(continuousData < 0) > 0
    tmpShiftVal = -1 * min(continuousData); 
else
    tmpShiftVal = 0;
end
continuousData = continuousData + tmpShiftVal;
max_conData = max(continuousData);

continuousData_scaled = continuousData / max_conData;

%
plot(verticalOffset + continuousData_scaled, 'LineWidth',llw, 'Color','k')

% hold all
tmpEventTime = find(eventData);
for iEvent = 1 : length(tmpEventTime)
    line([tmpEventTime(iEvent) tmpEventTime(iEvent)], verticalOffset + 1.2+[0 sl], ...
         'color','r', 'LineWidth',slw);
end

axis tight
axis off; box off


%% draft
% run v3_commonStuff.m

% %%
% signalParams.dt         = 1/signalParams.SF;               % length of time bin (in s)
% signalParams.t          = 0 : signalParams.dt : signalParams.signalLength - signalParams.dt;
% signalParams.nSample    = numel(signalParams.t);

% %%
% % this is useless
% % noiseLFP = ...
% %     globalDynamicsParams.whiteNoise_sigma ...
% %     * randn(signalParams.nCh, signalParams.nSample, signalParams.nTr);

% nOscComp = numel(globalDynamicsParams.oscFreq);
% % amps = globalDynamicsParams.oscComps;
% % oscWeights = amps / sum(amps);
% % oscFreqs = globalDynamicsParams.oscFreq;

% for iOscComp = 1 : nOscComp
%     cmplx_oscComps(iOscComp, :) =  ...
%         exp(1i * (2*pi* globalDynamicsParams.oscFreq(iOscComp) * signalParams.t));
% end

%%
% clear signalParams