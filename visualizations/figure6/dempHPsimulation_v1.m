%% prepare the data
ignit_gpla
load(fullfile(dirPath_gpla.data, 'dataHPsim', 'Session(full)_10_CA3CA1Network.mat'))

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

k = 1;
[~, sbh] = subplot2d(nR,nC, k, [3 4]);
displaceFigureStuff(sbh, [NaN k*upos+upos_cont NaN -upos])

plot(lfp(tr, itr, ich), 'k', 'linewidth', lw)
xlim([0 ns])

% plotScaleBar1(1:n(iDataset).LFPsamples, tlfp, '2 sec', NaN, 2 * Fs, NaN, ...
              % 0.05, NaN, 1.1, NaN, 1, 5.5); axis on


% title('raw LFP')
box off; axis off

k = 2;
[~, sbh] = subplot2d(nR,nC, k, [3 4]);
displaceFigureStuff(sbh, [NaN k*upos+upos_cont NaN -upos])

plot(lfpFilt_allSepFreq{iFreq}(ich, tr, itr), 'k', 'linewidth', lw)
xlim([0 ns])
% title('filtered LFP')
box off; axis off

k = 3;
[~, sbh_lastRowLfpDemo] = subplot2d(nR,nC, k, [3 4]);
displaceFigureStuff(sbh_lastRowLfpDemo, [NaN k*upos+.015+upos_cont NaN -upos])
tpfr = mean(full(spkTrains{itr}(:, tr)), 1);
plot(tpfr, 'k', 'linewidth', lw)
xlim([0 ns])
% title('population firing rate')
box off; axis off

