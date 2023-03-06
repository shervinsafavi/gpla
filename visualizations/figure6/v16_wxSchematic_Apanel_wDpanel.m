% the different with with v15, is some labels are different 

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
addpath /home/ssafavi/tools/utilities/matlab/functions/src/viz/

%% basics
clear all
ignit_gpla
ignit_datsrv(pds)


%%
% figure

%%
vc.f11.w = 17.4;
vc.f11.h = 23;

%%
fh = gcf;

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f11.w;
fh.Position(4)  = vc.f11.h;

%%
nR = 11-3;
nC = 4;

sdi = 3; % simulation demo index 

clf


%%
run('../vizConventions')

%% B
% demo traces
run dempHPsimulation_v4.m
pos_lastRowLfpDemo = get(sbh_lastRowLfpDemo, 'position');

%% C, E
% many of the stuff
run ca1_gpla_v7.m 
% toc

%% F
% spike-LFP phase diff
% mean filed 
% %% CA1-CA3
% tos = -.035 % tmp offset for next 2 subplots
% tos_lpp = -.065 % tmp offset for next 2 subplots (lpp : last polar plot)
tos_lpp = -.060 % tmp offset for next 2 subplots (lpp : last polar plot)
[~, sbh] = subplot2d(nR,nC, [sdi+4 sdi+5]-1, [1 2]);
pos_phaseDiff = get(sbh, 'Position'); 
yepd = .25*pos_phaseDiff(4); % yEnlarg_phaseDiff
% displaceFigureStuff(sbh, [NaN tos_lpp NaN NaN])
displaceFigureStuff(sbh, [NaN yepd+tos_lpp+.03 NaN -yepd])
% run ca1_neuralMass_v2.m
run ca1_neuralMass_EIphaseDiff_v0.m

% placeHolder
% xlabel('dummy')
% box on

%% G
% spike-LFP phase diff
% phase diff evolution
% tos_phEvo = -.03;
% tos_phEvo = -.05;
tos_phEvo = -.05;
% subplot2d(nR,nC, [sdi+3 sdi+4], [1 2])
[~, sbh] = subplot2d(nR,nC, [sdi+3]+2, [1 2]);
displaceFigureStuff(sbh, [NaN tos_phEvo NaN yepd]);

% run phaseEvolution_v2
% run phaseEvolution_v3.m
% we are not using phaseEvolution here anymore

% run ca1_neuralMass_v1.m
run ca1_neuralMass_spkLfpPhaseDiff_v0.m
ylabel({'Spike-LFP', 'Phase diff. [deg]'})

%% network schema
% *** to be placed manually
% [~, sbh] = subplot2d(nR,nC, [1 sdi], [1 2]);
% % displaceFigureStuff(sbh, [-.045 NaN .03 0.02])
% displaceFigureStuff(sbh, [-.045 .05 .03 0.02])
% % imshow('ca1ca3modelDiag.jpg');
% imshow('NeuronSchema_pdf.png');
% % imshow('NeuronSchema_pdf_v1.tiff');
% set(gca,'LooseInset', max(get(gca,'TightInset'), 0))


%% LFP part of ca1_gpla_v?
% need to be after "phase diff evolution"
run phaseEvolution_v3_onlyCompute.m
run ca1_gpla_lfpVec_v1.m
% run ca1_gpla_lfpVec_v1_wxImages.m

%%  CA1-CA3 joint analysis
[~, sbh] = subplot2d(nR,nC, [sdi+4 sdi+5], 2+[1 2])
displaceFigureStuff(sbh, [NaN tos_lpp NaN NaN])
run CA1CA3_gammaCoordination_v0.m

%% subplot labels
vc.f11.sblfz = 14;
% lr1y : left row 1 (y corrdinate)
vcrd.lr1y = 0.9641;
vcrd.lr2y = 0.69;
% vcrd.lr3y = 0.365;
vcrd.lr3y = 0.365;
% vcrd.lr4y = 0.22;
% vcrd.lr4y = 0.12;
vcrd.lr4y = 0.18;
% rr1y : right row 1 (y corrdinate)
vcrd.rr1y = vcrd.lr1y;
% vcrd.rr2y = 0.75; old
vcrd.rr2y = 0.77; % new
vcrd.rr3y = 0.605;
% vcrd.rr4y = vcrd.lr4y;
vcrd.rr4y = 0.22;


% vcrd.c1x = 0.01;
vcrd.c1x = 0;
% vcrd.c2x = 0.54;
vcrd.c2x = 0.51;

vcrd.all = [...
    vcrd.c1x vcrd.lr1y; ... % A
    vcrd.c2x vcrd.rr1y; ... % B
    vcrd.c2x vcrd.rr2y; ... % C
    vcrd.c1x vcrd.lr2y; ... % D
    vcrd.c2x vcrd.rr3y; ... % E
    vcrd.c1x vcrd.lr3y; ... % F
    vcrd.c1x vcrd.lr4y; ... % G
    vcrd.c2x vcrd.rr4y      % H
           ];

vcrd.labels = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H' };


for k = 1 : size(vcrd.all, 1)
    annotation('textbox',[vcrd.all(k,1) vcrd.all(k,2)  1 0.035], ...
               'String', vcrd.labels{k}, ...
               'FontSize', vc.f11.sblfz, ... 
               'LineWidth',1, ...
               'EdgeColor', 1*[1 1 1], ...
               'BackgroundColor', [1 1 1], ...
               'FitBoxToText','on')
end

%%
fsp = '/home/ssafavi/Nextcloud/research/nnr/reports_nnr/papers/p_ploscb2018_gpla/figures'
svflg = 1;
fn = 'fig11_v16_wxSchematic_Apanel_wDpanel_hpSimulation_extended';
spf(fsp, fn, svflg)


