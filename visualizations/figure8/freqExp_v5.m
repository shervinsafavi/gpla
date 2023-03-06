% in v4 lenend is improved visually

%%
% vc.ppdot = 6.5;
vc.ppdot = 3.5;
vc.c3r1f2.al = .35;
vc.c3r1f2.lw = 1;
vc.c3r1f2.fs = 5; % polar plots font size

% vc.f5.ppLW = 1.5; % moved to main vc file
vc.f5.ppFS = 4;

tos = -.05
txos = -.05; % displacment across x to creat space for colorbar

clear sbh
% clf
[~, sbh] = subplot2(nR,nC, iR, 1);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);

polarplot(angle(lfpVec.uate2010{ifreq}), ...
    abs(lfpVec.uate2010{ifreq}), 'k.', 'MarkerSize', vc.ppdot, ...
          'color',vc.lfpv);

% pax = gca;
% pax.FontSize = vc.c3r1f2.fs;
% pax.GridAlpha = vc.c3r1f2.al;
% pax.LineWidth = vc.c3r1f2.lw;
% thetaticks(0:45:315)


hold on
% subplot2(nR,nC, iR, 3);
polarplot(angle(spkVec.uate2010{ifreq}), ...
    abs(spkVec.uate2010{ifreq}), 'k.', 'MarkerSize', vc.ppdot, 'color',vc.spkv);

pax = gca;
pax.FontSize = vc.c3r1f2.fs;
pax.GridAlpha = vc.c3r1f2.al;
pax.LineWidth = vc.c3r1f2.lw;
thetaticks(0:45:315)

tr = get(gca, 'rlim');
trn = find_closestEvenNum1(tr(2));
rticks([0 trn/2 trn]);
rlim([0 trn])
pax.RAxisLocation = 67.5;

%%
% legend for spike and LFP vectors
if ifreq == 1
    [lh_pp icons] = legend('LFP vector', 'Spike vector');
    % displaceFigureStuff(lh_pp, [nan -.0615 nan nan]);
    displaceFigureStuff(lh_pp, [-.01 -.0523 nan nan]);
    % displaceFigureStuff(lh_pp, [-.25 .0621 nan nan]);
    set(lh_pp, 'box', 'off')
    set(lh_pp, 'fontsize', vc.f5.lfs)
    
    bof = +.1; % backward offset
    for k = 3 : 6
        tmpOrigPos = icons(k).XData;
        icons(k).XData = tmpOrigPos + bof; 
        icons(k).MarkerSize = 12;
    end

    for k = 1 : 2
        icons(k).FontSize = vc.f5.lfs;
    end

    % this is not needed
    % change the location of  text
    % for k = 1 : 3 % the first 3 are related text
    %     tmpOrigPos = icons(k).Position(1);
    %     icons(k).Position(1) = tmpOrigPos - bof; 
    % end
end

% axis off
% set(gca, 'fontsize', vc.f5.ppFS)
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%%
% subplot2(nR,nC, iR, 3);
clear sbh
[~, sbh] = subplot2(nR,nC, iR, 3);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);
% displaceFigureStuff(sbh, [txos tos NaN NaN]);

% image_mapArray(rad2deg(angle(lfpVec.uate2010{ifreq})), ...
%                utahMaps(iDataset).LFP.map);
tv = lfpVec.uate2010{ifreq};
image_mapArray(tv, ...
               utahMaps(iDataset).LFP.map);
for k = 0 : 10
   xline(k-.5);
   yline(k-.5);
end

% clim([-pi pi])
caxis([-180 180])
colormap hsv
axis on; box on; 
set(gca, 'ytick', {})
set(gca, 'xtick', {})

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%%
% subplot2(nR,nC, iR, 2);
[~, sbh] = subplot2(nR,nC, iR, 2);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);

theta = 0 : .01 : 2*pi;

svec = spkVec.uate2010{ifreq};
[mu, kappa] = getMacroPic(svec);
y = vmpdf1(theta, mu, kappa);

polarplot(theta, y * gPLV.uate2010(ifreq), ...
          'color',vc.spkv, ...
          'linewidth',vc.f5.ppLW);

svec = lfpVec.uate2010{ifreq};
[mu, kappa] = getMacroPic(svec);
y = vmpdf1(theta, mu, kappa);

hold on
polarplot(theta, y * gPLV.uate2010(ifreq), ...
          'color',vc.lfpv, ...
          'linewidth',vc.f5.ppLW);

pax = gca;
pax.FontSize = vc.c3r1f2.fs;
pax.GridAlpha = vc.c3r1f2.al;
pax.LineWidth = vc.c3r1f2.lw;
thetaticks(0:45:315)

% set(gca, 'fontsize', vc.f5.ppFS)
% axis off

tr = get(gca, 'rlim');
trn = find_closestEvenNum1(tr(2));
rticks([0 trn/2 trn]);
rlim([0 trn])
pax.RAxisLocation = 67.5;


set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%%
% we need this one to be last for the colorbar 

% subplot2(nR,nC, iR, 4);
[~, sbh] = subplot2(nR,nC, iR, 4);
% displaceFigureStuff(sbh, [NaN tos NaN NaN]);
displaceFigureStuff(sbh, [txos tos NaN NaN]);

% image_mapArray(rad2deg(angle(lfpVec.uate2010{ifreq})), ...
%                utahMaps(iDataset).LFP.map);
tv = spkVec.uate2010{ifreq};
image_mapArray(tv, ...
               utahMaps(iDataset).MU.map);
for k = 0 : 10
   xline(k-.5);
   yline(k-.5);
end

text(12, 5.5, [freqBandsLab{ifreq}, 'Hz'], ...
             'fontsize', 8.5)

% clim([-pi pi])
caxis([-180 180])
colormap hsv
axis on; box on; 
set(gca, 'ytick', {})
set(gca, 'xtick', {})

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
