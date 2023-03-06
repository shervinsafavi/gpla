%% basics
clear all
% ignit_gpla
pds = ignit()

%%
run('../vizConventions')

nR = 4;
nC = 5;

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

%% samples
subplot2(nR,nC, 1, 1+[1 2]);
subplot2(nR,nC, 1, 1+[3 4]);

%% gPLV across freqs
subplot2(nR,nC, 2, [1 2]);
% run gplaAcrossFreqs

%% phase difference
subplot2(nR,nC, 2, [3 4]);

%% freq 3-5
iR = 3;
ifreq = 1;
run freqExp

% subplot2(nR,nC, iR, 2);

%% freq 15 - 30 Hz
iR = 4;
ifreq = 3;
run freqExp

% subplot2(nR,nC, iR, 1);
% subplot2(nR,nC, iR, 2);

%%
fsp = '/home/ssafavi/Nextcloud/research/nnr/reports_nnr/papers/p_ploscb2018_gpla/figures'
svflg = 1;
fn = 'fig5_v2_UATE2010res';
spf(fsp, fn, svflg)
