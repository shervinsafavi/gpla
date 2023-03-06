vc.ppdot = 6.5;
vc.c3r1f2.al = .35;
vc.c3r1f2.lw = 1.5;
vc.c3r1f2.fs = 12;

subplot2(nR,nC, iR, 1);
polarplot(angle(lfpVec.uate2010{ifreq}), ...
    abs(lfpVec.uate2010{ifreq}), 'k.', 'MarkerSize', vc.ppdot, ...
          'color',vc.lfpv);

pax = gca;
pax.FontSize = vc.c3r1f2.fs;
pax.GridAlpha = vc.c3r1f2.al;
pax.LineWidth = vc.c3r1f2.lw;
thetaticks(0:45:315)


hold on
% subplot2(nR,nC, iR, 3);
polarplot(angle(spkVec.uate2010{ifreq}), ...
    abs(spkVec.uate2010{ifreq}), 'k.', 'MarkerSize', vc.ppdot, 'color',vc.spkv);

pax = gca;
pax.FontSize = vc.c3r1f2.fs;
pax.GridAlpha = vc.c3r1f2.al;
pax.LineWidth = vc.c3r1f2.lw;
thetaticks(0:45:315)


subplot2(nR,nC, iR, 2);
% image_mapArray(rad2deg(angle(lfpVec.uate2010{ifreq})), ...
%                utahMaps(iDataset).LFP.map);
image_mapArray(lfpVec.uate2010{ifreq}, ...
               utahMaps(iDataset).LFP.map);
axis on; box on; 
set(gca, 'ytick', {})
set(gca, 'xtick', {})

subplot2(nR,nC, iR, 3);
% image_mapArray(rad2deg(angle(lfpVec.uate2010{ifreq})), ...
%                utahMaps(iDataset).LFP.map);
image_mapArray(spkVec.uate2010{ifreq}, ...
               utahMaps(iDataset).MU.map);
axis on; box on; 
set(gca, 'ytick', {})
set(gca, 'xtick', {})

%%
subplot2(nR,nC, iR, 4);
theta = 0 : .01 : 2*pi;

svec = spkVec.uate2010{ifreq};
[mu, kappa] = getMacroPic(svec);
y = vmpdf1(theta, mu, kappa);

polarplot(theta, y * gPLV.uate2010(ifreq), ...
          'color',vc.spkv);

svec = lfpVec.uate2010{ifreq};
[mu, kappa] = getMacroPic(svec);
y = vmpdf1(theta, mu, kappa);

hold on
polarplot(theta, y * gPLV.uate2010(ifreq), ...
          'color',vc.lfpv);

