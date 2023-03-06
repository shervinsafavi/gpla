% difference with v5 is using different colormap for schematic in F


%%
% Each figure should be able to fit on a single 8.5 x 11 inch page. 
% Please do not send figure panels as individual files. We use three standard widths for figures: 
% 1 column, 85 mm; 
% 1.5 column, 114 mm; and 
% 2 column, 174 mm (the full width of the page). 
% Although your figure size may be reduced in the print journal, please keep these widths in mind. For Previews and other three-column formats, these widths are also applicable, though the width of a single column will be 55 mm.

clf

%%
% nR = 1+ 3*2 + 4;
nR = 1+ 3*1 + 4 + 1 - 1;
nC = 6;
nhc = nC / 2; % number half col

%
vc.f13.w = 17.4;
% vc.f13.w = 11.4;
vc.f13.h = vc.f13.w * nR / nC;

%
fh = gcf;

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f13.w;
fh.Position(4)  = vc.f13.h;

%%
hos = .03 * [-1 1]; % horizental offset
% vos_svs = -.07;     % vertical offset singular vectors
vos_svs = -.05;     % vertical offset singular vectors

vc.f13.fs_text = 6;

%%
pds = ignit();
% ignit_gpla
% ignit_datsrv

%% prep and data load
lfpIdx = [3,7,1,2];
lfpExt = {'Ve','exo','epsp','ipsp'};

% selBand =1:4;%1:length(globalDynamicsParams)
% selLFP = 3;

selBand = 1 : 4; %1:length(globalDynamicsParams)
nfreq = numel(selBand);
nModel = 2;

%%

% clf

eiNetImages = {'EInet_weakRecurrent.png', 'EInet_storngRecurrent.png'};

% vos_eiNets = .03 - 0.02;
vos_eiNets = .01;
hos_eiNets = .04;

% for kModel = 1 : nModel
%     if kModel == 1, coi = 0; 
%     elseif kModel == 2, coi = 3; end
%     [~, sbh] = subplot2d(nR,nC, 1, coi + nhc - 1)
%     displaceFigureStuff(sbh, [hos(kModel)-hos_eiNets vos_eiNets .1 .01])    
%     % imshow('EInet_weakRecurrent.png');

%     [A, map, transparency] = imread(eiNetImages{kModel});
%     image(A, 'AlphaData' , transparency);
%     % imshow (A, map, 'Parent' , gc)
%     % imshow (A)
%     axis off
%     pbaspect([1 .65 1])
%     % set(gca, 'Color' , [1 1 1] )
% end


%% example stuff for neural field model
% we run this line just to prepare for the line after
run demoSVecs_v2.m

%%
run demoNeuralFieldModel_v2.m

%%
% clf
netColors = cool(2);
% netColors = winter(2);
% vos_svs2 = -1.0700;
vos_svs2 = -0.0;
hos2 = 0.04;
hos3 = 0.06;
% hos2 = 0.06;
% hos3 = 0.08;
% trInd = [nR-1 nR]; % temp row index
% trInd = [nR-1 nR]; % temp row index
trInd = [4 5]; % temp row index
tmp = [1 2];
tcInd = [tmp; tmp+2; tmp+4];
run neuralFieldDemo_otherQuantities_v3.m
% subplot2d(nR, nC, [nR-1 nR], [1 2])
% subplot2d(nR, nC, [nR-1 nR], [1 2]+2)
% subplot2d(nR, nC, [nR-1 nR], [1 2]+4)

%% inset for schematic
% clf
inset([.81 .4 .2 .25])
pbaspect([.8 1 1])
run connKernelSchamatic_v4.m
% move the phase lag ti
displaceFigureStuff(th, [NaN -90 NaN])

%% circular colorbar - need to be below schematic
% clf
% inset([.8 .4 .2 .25]) from the other one
% inset([.835 .47 .095 .095]);
inset([.835 .36 .095 .095]);
axis square

% circ_collbar([0.1 1])
tickFontSize = vc.f13.fs_text;
% circ_collbar_wThTick1([0 1], tickFontSize)
tickLocs = ...
        [5.8 0; ...
         -.4 5.8; ...
         -6.7 0; ...
         -.0 -5.8]

circ_collbar_wThTick2([0 1], tickFontSize, tickLocs)
axis off
box off

%% subplot labels
vc.f13.sblfz = 14;
vcrd.r1y = 0.9691;
% vcrd.r2y = 0.5121;
vcrd.r2y = 0.565;
vcrd.r3y = 0.365;

vcrd.c1x = 0.00;
vcrd.c2x = 0.54
vcrd.c2x_r2 = .225;
vcrd.c2x_r4 = .68;
% vcrd.r3c2x = .33;
% vcrd.r3c3x = .61;

vcrd.all = [...
    vcrd.c1x vcrd.r1y; ... % A v
    vcrd.c2x vcrd.r1y; ... % B v
    vcrd.c1x vcrd.r2y; ... % C v
    vcrd.c2x vcrd.r2y; ... % E v 
    vcrd.c2x_r2 vcrd.r2y; ... % D v
    vcrd.c2x_r4 vcrd.r2y; ... % F v    
    vcrd.c1x vcrd.r3y; ... % G
    vcrd.c2x vcrd.r3y      % H
           ];

% vcrd.labels = {'A', 'B', 'C', 'D'};
vcrd.labels = {'A', 'B', 'C', 'E', 'D', 'F', 'G', 'H'};

for k = 1 : size(vcrd.all, 1)
    annotation('textbox',[vcrd.all(k,1) vcrd.all(k,2)  1 0.035], ...
               'String', vcrd.labels{k}, ...
               'FontSize', vc.f13.sblfz, ... 
               'LineWidth',1, ...
               'EdgeColor', 1*[1 1 1], ...
               'BackgroundColor', [1 1 1], ...
               'FitBoxToText','on')
end

% text(r2.x, 1.25, 'D', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize', vc.f13.sblfs)
% text(r2.x, r1.y, 'A', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize', vc.f13.sblfs)
% text(7, r1.y, 'B', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize', vc.f13.sblfs)
% text(14.25, r1.y, 'C', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize', vc.f13.sblfs)

%% cicualr colorbar
% % clf
% % inset([.8 .4 .2 .25]) from the other one
% % inset([.835 .47 .095 .095]);
% inset([.835 .35 .095 .095]);
% axis square

%% circular colorbar
% % circ_collbar([0.1 1])
% tickFontSize = vc.f13.fs_text;
% % circ_collbar_wThTick1([0 1], tickFontSize)
% tickLocs = ...
%         [5.8 0; ...
%          -.4 5.8; ...
%          -6.7 0; ...
%          -.0 -5.8]

% circ_collbar_wThTick2([0 1], tickFontSize, tickLocs)
% axis off
% box off



%%
% fsp = '/home/ssafavi/Nextcloud/research/nnr/reports_nnr/papers/p_ploscb2018_gpla/figures'
% svflg = 1;
% fn = 'fig13_v7_wxSchematic_ABpanel_neuralFieldModel';
% spf(fsp, fn, svflg)

%% draft
