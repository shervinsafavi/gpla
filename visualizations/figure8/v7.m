% in v7 we have only phase in the spatial maps
% and better allighment 

%%
% Each figure should be able to fit on a single 8.5 x 11 inch page. 
% Please do not send figure panels as individual files. We use three standard widths for figures: 
% 1 column, 85 mm; 
% 1.5 column, 114 mm; and 
% 2 column, 174 mm (the full width of the page). 
% Although your figure size may be reduced in the print journal, please keep these widths in mind. For Previews and other three-column formats, these widths are also applicable, though the width of a single column will be 55 mm.

%%
% later on should be added systematically
addpath /home/ssafavi/tools/utilities/matlab/toolboxes/utahArrays_preProcessing
addpath /home/ssafavi/tools/utilities/collections/matlab/functions/viz

%% basics
clear all
ignit_gpla
ignit_datsrv(pds)


%%
% figure

%%

vc.f5.w = 11.4;
vc.f5.h = 15;
vc.f5.fs_text = 5;
vc.f5.lfs = 5.5; % legend font size

%%
fh = gcf;

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f5.w;
fh.Position(4)  = vc.f5.h;


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
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%% LFP 
load([dirPath_gpla.data, filesep, 'lfpData_UA_20161019.mat'])
load([dirPath_gpla.data, filesep, 'metaData_UA_20161018.mat'])
% need it for chanel info

%% load data with analytical signals and whiten them
load(fullfile(dirPath_gpla.data, 'handies', ...
              'uate2010handyVars_wAnalyticLfp.mat'))

% too slow to load from remote server
% load(fullfile(pds.datsrv.hnd, ...
%               'uate2010handyVars_wAnalyticLfp.mat'))

% iDataset = 1;

%% LFP 
% *** should be moved toward the end for memory issues

% clf %
subplot2(nR,nC, 1, [2 4]);

% subplot(2,2, [1 2])

% nlfpt = 7;

% loi = [69 7 69 72 69 74 88 89 3 4];

% loi = [17 94 18 90 77 81 56 27 34 88]
loi = [17 94 18 90 77 81 54 25 70 41];

% loi = diag(utahMaps(iDataset).LFP.map);
% loi(2) = 47;
% loi(end) = 47;

nlfpt = numel(loi);
% loi = 68 * ones(1,nlfpt)

[elecsRow, elecsCol] = channelLocFinder2(loi, utahMaps(iDataset).LFP.map);


k = 1;
tlfp = lfp(iDataset).v(loi(k), :, iTr) + abs(min(lfp(iDataset).v(loi(k), :, iTr)));

plotScaleBar1(1:n(iDataset).LFPsamples, tlfp, '2 sec', NaN, 2 * Fs, NaN, ...
              0.05, NaN, 1.1, NaN, 1, 5.5); axis on

tlfp_old = tlfp;
hold all
for k = 1 : nlfpt
    tlfp = lfp(iDataset).v(loi(k), :, iTr) + abs(min(lfp(iDataset).v(loi(k), :, iTr)));
    % lfp2plot = tlfp_old + 1.1 * max(tlfp_old);
    lfp2plot = tlfp + max(tlfp_old);
    % tlfp = lfp(iDataset).v(loi(k), :, iTr);
    % lfp5plot = tlfp;
    % plot(1.1*max(lfp(iDataset).v(k-1, :, iTr)) + tlfp, 'k');
    % lfp2plot = 3*max(tlfp_old) + tlfp;
    plot(lfp2plot, 'k');
    tlfp_old = lfp2plot;
    % hold on
    text(-270, mean(lfp2plot), num2str(loi(k)), ...
         'HorizontalAlignment','center', 'fontsize',vc.f5.fs_text);
end

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

axis tight
box off
axis off

%% array map
% subplot 223
tos = .011;
[~, sbh] = subplot2(nR,nC, 2, 1);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);

% imagescnan(NaN(10));
% axis on
% box 
for k = 0 : 10
   xline(k);
   yline(k);
end
xlim([0 10])
ylim([0 10])


hold all
for k = 1 : nlfpt
    text(elecsCol(k)-.5, elecsRow(k)-.5, num2str(loi(k)), ...
         'HorizontalAlignment','center', 'fontsize',vc.f5.fs_text);
end
% grid on
axis square
axis off

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

% utahMaps(iDataset).LFP.map %

%% sample spike and LFP 
% we rather use LFP in different band
% right after LFP, to have LFP in memory
% subplot2(nR,nC, nresl, [1 2]);
% run demo_spkLfpRaw_v0.m

%% spike 
% clf
[~, sbh] = subplot2(nR,nC, 2, [2 4]);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);
%
% clf
plot_eventRaster(full(sparseSPT{iTr}), {'color', 'k'});
% plotScaleBar1([0 n(iDataset).LFPsamples], [1 n(iDataset).MU], '2 sec', NaN, 2 * Fs, NaN, ...
%               -0.6, NaN, -.8, NaN, 1, 5.5); axis on
% hold on
plotScaleBar1([0 n(iDataset).LFPsamples], [1 n(iDataset).MU], '2 sec', NaN, 2 * Fs, NaN, ...
              -1.1, NaN, -.1, NaN, 1, 5.5); axis on


axis off

%% sample spike and filtered LFP low freq

% clf
tos = .014;
tr = 800 : 1600;
[~, sbh] = subplot2(nR,nC, nresl, [1 2]);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);
ifreq = 1;
run demo_spkLfpFiltered_v1.m
plotScaleBar1(1:numel(tr), [0 4], '500 ms', NaN, Fs / 2, NaN, ...
              0.05, NaN, .1, NaN, 1, 5.5); axis on

axis off

% sample spike and filtered LFP beta 
tr = 500 : 700;
[~, sbh] = subplot2(nR,nC, nresl, 2+[1 2]);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);
ifreq = 3;
run demo_spkLfpFiltered_v1.m
plotScaleBar1(1:numel(tr), [0 4], '100 ms', NaN, Fs / 10, NaN, ...
              0.05, NaN, .1, NaN, 1, 5.5); axis on

axis off

% subplot2(nR,nC, nresl, [1 2]);
% subplot2(nR,nC, nresl+1, [1 2]);
% subplot2(nR,nC, nresl, 2+[1 2]);
% subplot2(nR,nC, nresl+1,2+ [1 2]);
% subplot2(nR,nC, 1, 1+[3 4]);

%% gPLV across freqs

[~, sbh] = subplot2(nR,nC, nred+1, [1 2]);
run gplaAcrossFreqs_v1.m
% the intial part should go into a seperate scripts (for loading
% the data)
% run gplaAcrossFreqs

%% phase gradient
ifreq = 3;
% delete(sbh)
[~, sbh] = subplot2(nR,nC, nred+1, [3 4]);
% displaceFigureStuff(sbh, [.01 .03 -.01 NaN])
tos = .03; % temp offset
displaceFigureStuff(sbh, [tos -.03 -tos NaN]);
run calcLfpPhaseGradient_v0.m

%% freq 15 - 30 Hz
% if we run the 3-5 Hz first it will overlap
iR = nred+3;
ifreq = 3;
run freqExp_v3.m
pos_immageArray2 = get(sbh, 'position');


% freq 3-5
iR = nred+2;
ifreq = 1;
run freqExp_v3.m
pos_immageArray1 = get(sbh, 'position');

% subplot2(nR,nC, iR, 2);


% subplot2(nR,nC, iR, 1);
% subplot2(nR,nC, iR, 2);
% colorbar

%%
cbh = colorbar;
set(cbh, 'position', [...
    pos_immageArray2(1) + pos_immageArray2(3) + .005 ...
    pos_immageArray2(2) ...
    .022 ...
    2*pos_immageArray1(4) + (pos_immageArray1(2) - (pos_immageArray2(2) + pos_immageArray2(4))) ...
                   ])

title(cbh, 'Phase [angle]')
set(cbh, 'fontsize', vc.f5.lfs)

%%
fsp = '/home/ssafavi/Nextcloud/research/nnr/reports_nnr/papers/p_ploscb2018_gpla/figures'
svflg = 1;
fn = 'fig5_v7_UATE2010res';
spf(fsp, fn, svflg)
