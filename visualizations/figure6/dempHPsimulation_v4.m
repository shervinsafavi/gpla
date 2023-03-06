% in v3 we have extra labels

%% prepare the data
% ignit
pds = ignit();
load(fullfile(pds.ldat, 'dataHPsim', 'Session(full)_10_CA3CA1Network.mat'))
% load(fullfile(dirPath_gpla.data, 'dataHPsim', 'Session(full)_10_CA3CA1Network.mat'))

%%
spk = Sig.spikes;
spkT = Sig.timespk;
lfp = Sig.datca1;
% lfp = cat(3,Sig.datca1,Sig.datca3);
lfpT = Sig.time;
lfpT = lfpT+.5;% align reference with spikes

% resample spikes at lfp resolution
spkSparse = zeros(length(lfpT),size(lfp,2),size(spk,3));
for kTime = 1:length(spkT)
    spkSparse(floor(spkT(kTime)/Sig.dx)+1,:,:)=spkSparse(floor(spkT(kTime)/Sig.dx)+1,:,:)+spk(kTime,:,:);
end

spkSparse = logical (spkSparse);
spkTrains = {};
for kTr = 1:size(spkSparse,2)
    spkTrains{kTr} = sparse(squeeze(spkSparse(:,kTr,:))');
end

[signalParams.nSample, signalParams.nTr, signalParams.nCh] = size(lfp);

clear lfpPh_allSpeFreq

Bands = {[1,5],[5,20],[20,40],[40,80],[80,180]};

nFreq = numel(Bands);
for iFreq = 1 : nFreq 
   freqBandsLab{iFreq} = mat2str(Bands{iFreq});
end
freqBandsLab = {...
    '1-5';...
    '5-20'; ...
    '20-40'; ...
    '40-80'; ...
    '80-180' ...
    };

Fs = 1 / Sig.dx;
n = signalParams;

%% filter
filterOrder = 3;


% for iFreq = 1:length(Bands)
for iFreq = 5
    [b, a] = butter(filterOrder, Bands{iFreq}*Sig.dx*2, 'bandpass');
    for iTr = 1 : signalParams.nTr
        for iCh = 1 : signalParams.nCh
            tmpLfp = lfp(:, iTr,iCh);
            tmpLfp_filtered = filtfilt(b, a, tmpLfp);
            lfpFilt_allSepFreq{iFreq}(iCh,:, iTr) = tmpLfp_filtered;
            tmpAnlcLfp = hilbert(tmpLfp_filtered);
            lfpAnalSig_allSepFreq{iFreq}(iCh,:, iTr) = tmpAnlcLfp;
        end
    end
end

freqBands = Bands;

%%
ich = 1;
itr = 1;
% tr = 45 : 140;
tr = 1 : 340;
ns = numel(tr);

lw = 1.4;

upos = .04; % upward offset
upos_cont = .03; % upward offset
% upos = .06; % upward offset

% clf

exlim = 25; % extra for xlim
xenlarge = .05; % enlarge  xaxis
xslos = .2;
yslos = .02;
xstos = .25;
ystos = .01;
stoslw = 1;
slosfs = 6.5;

titos = .3;

k = 1;
[~, sbh] = subplot2d(nR,nC, k, [3 4]);
% displaceFigureStuff(sbh, [NaN k*upos+upos_cont NaN -upos])
displaceFigureStuff(sbh, [NaN k*upos+upos_cont xenlarge -upos])
%

tts = lfp(tr, itr, ich); % temp time series
plot(tts, 'k', 'linewidth', lw)
xlim([0 ns+exlim])

% plotScaleBar1(tr,tts, '50 ms', 'random y', 50, 10e-3, .2, .02, .25, .02, 1, 5.5); axis on
plotScaleBar1(tr,tts, '50 ms', '0.1 mV', 50, 10e-3, xslos, yslos, xstos, ystos, stoslw, slosfs); axis on
% xlim([0 500])

% title('raw LFP')
th = title('Raw LFP');
set(th, 'fontweight','normal')
set(th, 'fontsize', slosfs)
tpos = get(th, 'position');
tpos(1) = 0 + tpos(1)*titos;
set(th, 'position', tpos);
box off; axis off
% box off; axis off


%
k = 2;
[~, sbh] = subplot2d(nR,nC, k, [3 4]);
% displaceFigureStuff(sbh, [NaN k*upos+upos_cont NaN -upos])
displaceFigureStuff(sbh, [NaN k*upos+upos_cont xenlarge -upos])

plot(lfpFilt_allSepFreq{iFreq}(ich, tr, itr), 'k', 'linewidth', lw)
plotScaleBar1(tr,tts, '50 ms', '0.05 mV', 50, 5e-3, xslos/1.6, yslos, xstos/2.5, ystos, stoslw, slosfs); axis on
% xlim([0 ns])
xlim([0 ns+exlim])
% title('filtered LFP')
th = title('20-40 Hz filtered LFP');
set(th, 'fontweight','normal')
set(th, 'fontsize', slosfs)
tpos = get(th, 'position');
% tpos(1) = 0;
tpos(1) = 0 + tpos(1)*titos;
set(th, 'position', tpos);
box off; axis off
% box off; axis off

k = 3;
[~, sbh_lastRowLfpDemo] = subplot2d(nR,nC, k, [3 4]);
% displaceFigureStuff(sbh_lastRowLfpDemo, [NaN k*upos+.015+upos_cont NaN -upos])
displaceFigureStuff(sbh_lastRowLfpDemo, [NaN k*upos+.015+upos_cont xenlarge -upos])
tpfr = sum(full(spkTrains{itr}(:, tr)), 1);
plot(tpfr, 'k', 'linewidth', lw)
% plotScaleBar1(tr,tts, '50 ms', '0.2 mV', 50, 20e-3, xslos*1.2, yslos, xstos*1.7, ystos, stoslw, slosfs); 
% plotScaleBar1(tr,tts, '50 ms', '0.2 mV', 50, 20e-3, xslos*2, yslos, xstos*1.7, ystos, stoslw, slosfs); 
%%
plotScaleBar1(tr,tts, '50 ms', '5 spikes', 50, 5, xslos*500, yslos, xstos*500.7, ystos, stoslw, slosfs); 
% axis on
% xlim([0 ns])
xlim([0 ns+exlim])
th = title('Population firing rate');
set(th, 'fontweight','normal')
set(th, 'fontsize', slosfs)
tpos = get(th, 'position');
tpos(1) = 0 + tpos(1)*titos;
% tpos(1) = 0;
set(th, 'position', tpos);
box off; axis off

