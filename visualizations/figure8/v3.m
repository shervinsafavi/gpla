%% basics
clear all
% ignit_gpla
pds = ignit()

%%
run('../vizConventions')

nR = 7;
nC = 4;
nred = 4; % nRow example data

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

%% LFP 
subplot2(nR,nC, 1, [2 4]);


%%
subplot2(nR,nC, 2, [2 4]);
nresl = 3; % nRow example spike LFP 
subplot2(nR,nC, nresl, [1 2]);
subplot2(nR,nC, nresl+1, [1 2]);
subplot2(nR,nC, nresl, 2+[1 2]);
subplot2(nR,nC, nresl+1,2+ [1 2]);
% subplot2(nR,nC, 1, 1+[3 4]);

%% gPLV across freqs
subplot2(nR,nC, nred+1, [1 2]);
% run gplaAcrossFreqs

%% phase difference
subplot2(nR,nC, nred+1, [3 4]);

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
fn = 'fig5_v3_UATE2010res';
spf(fsp, fn, svflg)
