%% Figure 3: Statitical analysis of GPLA


%% initiatation 
% (e.g add necessary packages and functions)
clear all
clf
pds = ignit()

%% prepare the figure

% visualization conventions
vc = get_vizConventions();

%
fh = gcf;

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f10.w;
fh.Position(4)  = vc.f10.h;


%% assign some parameters 
run v6_commonStuff.m

%% subplot labels
vcrd.r1y = 0.9541;
vcrd.r2y = 0.675;
vcrd.r3y = 0.328;

vcrd.c1x_r1 = 0.00;
vcrd.c2x_r1 = 0.3;
vcrd.c3x_r1 = 0.682;

vcrd.c1x_r2 = vcrd.c1x_r1;
vcrd.c2x_r2 = 0.5;

vcrd.c1x_r3 = vcrd.c1x_r1;
vcrd.c2x_r3 = 0.27;
vcrd.c3x_r3 = 0.73;


vcrd.c2x = .286;
vcrd.c3x = .50;
vcrd.c4x = .705;

vcrd.all = [...
    vcrd.c1x_r1 vcrd.r1y; ... % A
    vcrd.c2x_r1 vcrd.r1y; ... % B
    vcrd.c3x_r1 vcrd.r1y; ... % C
    vcrd.c1x_r2 vcrd.r2y; ... % D
    vcrd.c2x_r2 vcrd.r2y; ... % E
    vcrd.c1x_r3 vcrd.r3y; ... % F
    vcrd.c2x_r3 vcrd.r3y; ... % G
    vcrd.c3x_r3 vcrd.r3y      % H
           ];

vcrd.labels = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};


for k = 1 : size(vcrd.all, 1)
    annotation('textbox',[vcrd.all(k,1) vcrd.all(k,2)  1 0.035], ...
               'String', vcrd.labels{k}, ...
               'FontSize', vc.f10.sblfz, ... 
               'LineWidth',1, ...
               'EdgeColor', 1*[1 1 1], ...
               'BackgroundColor', [1 1 1], ...
               'FitBoxToText','on')
end

%% LFP synthetiz
r1_hos = -.005
[~, sbh] = subplot2d(nR,nC, 1, 1);
displaceFigureStuff(sbh, [r1_hos .025 NaN NaN])
run v0_demoLfpSynthz.m

tr = 1 : 500;

oscCompId = nOscComp : -1 : 1;

for iOscComp = 1 : nOscComp
    tosc = real(3*iOscComp + cmplx_oscComps(iOscComp, tr));
    plot(tosc, 'k', 'linewidth', 2);
    txt = ['$$O_', num2str(oscCompId(iOscComp)), '(t)$$'];   
    text(-200, mean(tosc), txt, 'Interpreter', 'latex')
    hold on
end

th = title('Oscillatory components')

set(th, 'fontweight', 'normal')
displaceFigureStuff(th, [NaN 1 NaN])

axis tight

xline(403, 'r:', 'linewidth',2);

txt = '$$LFP(t) = \sum_{k = 1}^{5} w_k O_k(t)$$';   

axis off 
box off

%% spike LFP coupling demo
% probably need to be ran later
% run v0_demoLfpSpk.m

lr_vos = -.03;
r2_vos = .015;

%% detection as a function of coupling and dimention
caseName = 'correctPosInves'; 
caseNames{1} = caseName;

% run v2_couplingDetection
% fn = fullfile(storagePath, 'summaryStats_couplingDetection.mat');
% save(fn, 'summaryStat', '-v7.3');

load(fullfile('.', 'dat', 'summaryStats_couplingDetection'))

[~, sbh] = subplot2d(nR,nC, 3, 1);
displaceFigureStuff(sbh, [-.03 lr_vos NaN NaN])

tmpImageData = (sum(summaryStat.correctPosInves, 3))' / nRel;

imagesc(unitNums, (1:nCouplingStrength), tmpImageData * 100);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
xlabel('Num. of units/LFPs')
ylabel('Coupling strength')

ylim([2 10])
yticks((2:nCouplingStrength))
yticklabels(couplingStrengths(2:end))

colormap(gca, 'copper')
c = colorbar; 
title(c, 'Detection [%]')
set(gca,'Ydir','Normal')

box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [vc.f10.imts 0.035];

set(gca, 'fontsize', vc.f10.gfs)

%% check the false positive
caseName = 'falsePosInves'; 
caseNames{2} = caseName;

% run v2_falsePositive
% fn = fullfile(storagePath, 'summaryStats_falsePosInves.mat');
% save(fn, 'summaryStat', '-v7.3');

load(fullfile('.', 'dat', 'summaryStats_falsePosInves.mat'))

barData = sum(summaryStat.falsePosInves, 2);
[~, sbh] = subplot2d(nR,nC, 3, [2 3]);
displaceFigureStuff(sbh, [NaN lr_vos NaN NaN])

bar(unitNums, barData, 'k')
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
tfh = yline(5, 'g');
ylim([0 10])
yticks([0 5 10])
xlim([unitNums(1)-50 unitNums(end)+50])
xtickangle(45)
xlabel('Number of spiking units/LFP channels')
ylabel('Detection [%]')
grid on

set(gca, 'fontsize', vc.f10.gfs)

box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.035];


lh = legend(tfh, 'Significance threshold');
set(lh, 'fontsize', vc.f10.lfs)

%% multi-population detection
caseName = 'multPopDetection'; 
caseNames{3} = caseName;

% run v2_multPopDetection
% fn = fullfile(storagePath, 'summaryStats_multPopDetection.mat');
% save(fn, 'summaryStat', '-v7.3');

run v3_commonStuff
load(fullfile('.', 'dat', 'summaryStats_multPopDetection'))

clear tmpImageData
tmpImageData = NaN(nCouplingStrength, nPopNum);
for ipn = 1 : nPopNum
    for ic = 2 : nCouplingStrength
        tv = squeeze(summaryStat.multPopDetection(ipn, ic, :));
        tmpImageData(ic, ipn) = mean(((tv - ipn)/ipn) .^ 2);
    end
end

ic = 1
for ipn = 1 : nPopNum
    tv = squeeze(summaryStat.multPopDetection(ipn, ic, :));
    % tmpImageData(ic, ipn) = mean((tv - 0) .^ 2);
    tmpImageData(ic, ipn) = NaN;
end

[~, sbh] = subplot2d(nR,nC, 3, 4);
displaceFigureStuff(sbh, [.08 lr_vos NaN NaN])
imagesc(popNums, (1:nCouplingStrength), tmpImageData);
ylim([2 10])

yticks((2:nCouplingStrength))
yticklabels(couplingStrengths(2:end))

xlabel('Num. of population')
ylabel('Coupling strength')

colormap(gca, 'autumn')
c = colorbar; title(c, 'MSE'); 
set(gca,'Ydir','Normal')

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
set(gca, 'fontsize', vc.f10.gfs)

box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [vc.f10.imts 0.035];

% run this two at the end, as some of the parameters are
% overwritten e.g. nRel

%% MP distribution - with coupling
% tic; run v5_MPdistribution_wCoupling; toc
[~, sbh] = subplot2(nR,nC, 2, [3 4]);
displaceFigureStuff(sbh, [-r1_hos r2_vos NaN NaN])
% displaceFigureStuff(sbh, [-r1_hos .05 NaN NaN])
run v8_MPdistribution_wCoupling_plot.m
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
set(gca, 'fontsize', vc.f10.gfs)

% set the same lim for the other subplot
% subplot2(nR,nC, 2, [3 4]);
% currentYlim = get(gca, 'ylim')
% tyl = find_closestEvenNum1(currentYlim(2));
tyl = .16;
currentYlim = [0 tyl];
ylim([0 tyl])
yticks([0 tyl/2 tyl])
xlim([0 5])

%% MP distribution - no coupling
% tic; run v3_MPdistribution_noCoupling.m; toc
[~, sbh] = subplot2(nR,nC, 2, [1 2]);
displaceFigureStuff(sbh, [r1_hos r2_vos NaN NaN])
run v6_MPdistribution_noCoupling_plot.m
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
set(gca, 'fontsize', vc.f10.gfs)

ylim(currentYlim)
yticks([0 tyl/2 tyl])
vline(lambda, 'b')

%% REPEATED: MP distribution - with coupling
% some issues with variables
% tic; run v5_MPdistribution_wCoupling; toc
[~, sbh] = subplot2(nR,nC, 2, [3 4]);
displaceFigureStuff(sbh, [-r1_hos r2_vos NaN NaN])
run v8_MPdistribution_wCoupling_plot.m
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
set(gca, 'fontsize', vc.f10.gfs)

[vlh, txth] = vline(lambda, 'b', 'Theoretical threshold');
displaceFigureStuff(txth, [NaN .055 NaN])
set(txth, 'fontsize', vc.f10.lfs+.5)

% set the same lim for the other subplot
% subplot2(nR,nC, 2, [3 4]);
% currentYlim = get(gca, 'ylim')
% tyl = find_closestEvenNum1(currentYlim(2));
tyl = .16;
currentYlim = [0 tyl];
ylim([0 tyl])
yticks([0 tyl/2 tyl])
xlim([0 5])

%% example coupling matrix 
[~, sbh] = subplot2d(nR,nC, 1, 4);
[lfpLikeSig, spikeTrains] = ...
    smlt_sustLockedSpkLfpPairs_multFreqCoupling2(...
        globalDynamicsParams, spikeTrainParams, couplingParams, ...
        signalParams);

% filter LFP 
lfpPh = tpp_filt_recenter(lfpLikeSig, freqCenter, halfFilterWidth, ...
                                 signalParams.SF, filterOrder, nIteCent);

% GPLA

[~,~,~,~, ~, tmpSvdOut] = ...
    tngpla(spikeTrains, lfpPh, [], [], [], [] , [], ...
           statTestInfo, 'all', sameElecCheckInfo, ...
           'nSpk', 0, 0);

imagesc(abs(tmpSvdOut.couplingMatrix)); 
colormap(gca, 'copper')
c = colorbar; title(c, 'Coupling strength'); 
displaceFigureStuff(c, [0.08 NaN NaN NaN])
ylabel('LFP ID')
xlabel('Spiking unit ID')
axis square
set(gca,'Ydir','Normal')
set(gca, 'fontsize', vc.f10.gfs)
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%% spike LFP coupling demo
[~, sbh] = subplot2(nR,nC, 1, [2 3]);
% displaceFigureStuff(sbh, [r1_hos NaN NaN NaN])
displaceFigureStuff(sbh, [-.016 NaN NaN NaN])

run v0_demoLfpSpk.m


%%
box off
ax = gca;
ax.TickDir = 'out';
ax.TickLength = [vc.f10.imts 0.035];
