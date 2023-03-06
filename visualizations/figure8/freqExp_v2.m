vc.ppdot = 6.5;
vc.c3r1f2.al = .35;
vc.c3r1f2.lw = 1;
vc.c3r1f2.fs = 5; % polar plots font size

vc.f5.ppLW = 1.5;
vc.f5.ppFS = 4;

tos = -.05

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

% axis off
% set(gca, 'fontsize', vc.f5.ppFS)
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%
% subplot2(nR,nC, iR, 3);
clear sbh
[~, sbh] = subplot2(nR,nC, iR, 3);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);


% image_mapArray(rad2deg(angle(lfpVec.uate2010{ifreq})), ...
%                utahMaps(iDataset).LFP.map);
image_mapArray(lfpVec.uate2010{ifreq}, ...
               utahMaps(iDataset).LFP.map);
axis on; box on; 
set(gca, 'ytick', {})
set(gca, 'xtick', {})

set(gca,'LooseInset', max(get(gca,'TightInset'), 0))


% subplot2(nR,nC, iR, 4);
[~, sbh] = subplot2(nR,nC, iR, 4);
displaceFigureStuff(sbh, [NaN tos NaN NaN]);


% image_mapArray(rad2deg(angle(lfpVec.uate2010{ifreq})), ...
%                utahMaps(iDataset).LFP.map);
image_mapArray(spkVec.uate2010{ifreq}, ...
               utahMaps(iDataset).MU.map);
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

