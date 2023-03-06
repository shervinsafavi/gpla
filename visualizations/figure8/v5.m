%% basics
clear all
% ignit_gpla
pds = ignit()

%%
% load(fullfile('dat','gplaAcrossFreqs_v0'))
% iDataset = 1;

%%
run('../vizConventions')

nR = 6;
nC = 4;
nred = 3; % nRow example data
nresl = 3; % nRow example spike LFP 

iTr = 1;

clf 

%%
freqBands = [...
    3 5;...
    5 15; ...
    15 30; ...
    30 50; ...
    50 70; ...
    70 105; ...
    105 125; ...
    125 140, ...
            ];
freqBandsLab = {...
    '3-5';...
    '5-15'; ...
    '15-30'; ...
    '30-50'; ...
    '50-70'; ...
    '70-105'; ...
    '105-125'; ...
    '125-140', ...
               };

% iDataset = 1

%%
subplot2(nR,nC, 1, 1);
imshow('UtahRecSite.png')

%%
subplot2(nR,nC, 2, 1);
imagescnan(NaN(10))
grid on
axis square 


%% LFP 
load([dirPath_gpla.data, filesep, 'lfpData_UA_20161019.mat'])

%%
subplot2(nR,nC, 1, [2 4]);
% clf
% nlfpt = 7;
loi = [69 7 69 72 69 74 88 89 3 4];
nlfpt = numel(loi)
% loi = 68 * ones(1,nlfpt)

k = 1;
tlfp = lfp(iDataset).v(loi(k), :, iTr) + abs(min(lfp(iDataset).v(loi(k), :, iTr)));
tlfp_old = tlfp;
for k = 2 : nlfpt
    tlfp = lfp(iDataset).v(loi(k), :, iTr) + abs(min(lfp(iDataset).v(loi(k), :, iTr)));
    % lfp2plot = tlfp_old + 1.1 * max(tlfp_old);
    lfp2plot = tlfp + max(tlfp_old);
    % tlfp = lfp(iDataset).v(loi(k), :, iTr);
    % lfp5plot = tlfp;
    % plot(1.1*max(lfp(iDataset).v(k-1, :, iTr)) + tlfp, 'k');
    % lfp2plot = 3*max(tlfp_old) + tlfp;
    plot(lfp2plot, 'k');
    tlfp_old = lfp2plot;
    hold on
end

axis tight
box off
axis off

%% sample spike and LFP 
% we rather use LFP in different band
% right after LFP, to have LFP in memory
% subplot2(nR,nC, nresl, [1 2]);
% run demo_spkLfpRaw_v0.m

%% spike 
subplot2(nR,nC, 2, [2 4]);
plot_eventRaster(full(sparseSPT{iTr}), {'color', 'k'});

%% sample spike and filtered LFP low freq
tr = 800 : 1600;
subplot2(nR,nC, nresl, [1 2]);
ifreq = 1;
run demo_spkLfpFiltered_v0

%% sample spike and filtered LFP beta 
% tr = 500 : 700;
tr = 800 : 700;
subplot2(nR,nC, nresl, 2+[1 2]);
ifreq = 3;
run demo_spkLfpFiltered_v0

% subplot2(nR,nC, nresl, [1 2]);
% subplot2(nR,nC, nresl+1, [1 2]);
% subplot2(nR,nC, nresl, 2+[1 2]);
% subplot2(nR,nC, nresl+1,2+ [1 2]);
% subplot2(nR,nC, 1, 1+[3 4]);

%% gPLV across freqs
subplot2(nR,nC, nred+1, [1 2]);
run gplaAcrossFreqs_v0
% the intial part should go into a seperate scripts (for loading
% the data)
% run gplaAcrossFreqs

%% phase gradient
subplot2(nR,nC, nred+1, [3 4]);
run calcLfpPhaseGradient_v0

%% freq 3-5
iR = nred+2;
ifreq = 1;
run freqExp_v1

% subplot2(nR,nC, iR, 2);

%% freq 15 - 30 Hz
iR = nred+3;
ifreq = 3;
run freqExp_v1

% subplot2(nR,nC, iR, 1);
% subplot2(nR,nC, iR, 2);

%%
fsp = '/home/ssafavi/Nextcloud/research/nnr/reports_nnr/papers/p_ploscb2018_gpla/figures'
svflg = 1;
fn = 'fig5_v5_UATE2010res';
spf(fsp, fn, svflg)
