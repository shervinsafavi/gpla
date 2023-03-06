%% Figure 6: GPLA on HP-SWR simulations

%% initiatation 
% (e.g add necessary packages and functions)
clear all
clf
pds = ignit();


%% prepare the figure

% visualization conventions
vc = get_vizConventions();

fh = gcf;

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f11.w;
fh.Position(4)  = vc.f11.h;

% number of rows and columns in the subplot
nR = 11-3;
nC = 4;

sdi = 3; % simulation demo index

%%
vc = get_vizConventions();

%% B
% demo traces
run dempHPsimulation_v4.m
pos_lastRowLfpDemo = get(sbh_lastRowLfpDemo, 'position');

%% C, E
run ca1_gpla_v7.m

%% F
% spike-LFP phase differences
% mean filed 
tos_lpp = -.060 % tmp offset for next 2 subplots (lpp : last polar plot)
[~, sbh] = subplot2d(nR,nC, [sdi+4 sdi+5]-1, [1 2]);
pos_phaseDiff = get(sbh, 'Position'); 
yepd = .25*pos_phaseDiff(4); % yEnlarg_phaseDiff
displaceFigureStuff(sbh, [NaN yepd+tos_lpp+.03 NaN -yepd])
% run ca1_neuralMass_v2.m
run ca1_neuralMass_EIphaseDiff_v0.m

%% G
% spike-LFP phase diff
% phase differece evolution
tos_phEvo = -.05;
[~, sbh] = subplot2d(nR,nC, [sdi+3]+2, [1 2]);
displaceFigureStuff(sbh, [NaN tos_phEvo NaN yepd]);

% run phaseEvolution_v2
% run phaseEvolution_v3.m

% run ca1_neuralMass_v1.m
run ca1_neuralMass_spkLfpPhaseDiff_v0.m
ylabel({'Spike-LFP', 'Phase diff. [deg]'})

%% network schema
[~, sbh] = subplot2d(nR,nC, [1 sdi], [1 2]);
displaceFigureStuff(sbh, [-.045 .05 .03 0.02])
imshow('NeuronSchema_pdf.png');
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))


%% LFP part of ca1_gpla_v?
% need to be after "phase diff evolution"
run phaseEvolution_v3_onlyCompute.m
run ca1_gpla_lfpVec_v1.m

%%  CA1-CA3 joint analysis
[~, sbh] = subplot2d(nR,nC, [sdi+4 sdi+5], 2+[1 2])
displaceFigureStuff(sbh, [NaN tos_lpp NaN NaN])
run CA1CA3_gammaCoordination_v0.m

%% subplot labels
vc.f11.sblfz = 14;
vcrd.lr1y = 0.9641;
vcrd.lr2y = 0.69;
vcrd.lr3y = 0.365;
vcrd.lr4y = 0.18;
vcrd.rr1y = vcrd.lr1y;
vcrd.rr2y = 0.77; % new
vcrd.rr3y = 0.605;
vcrd.rr4y = 0.22;


vcrd.c1x = 0;
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

