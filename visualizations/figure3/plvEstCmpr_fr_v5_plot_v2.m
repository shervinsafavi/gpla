% in v2, we plot less firing rates

% clf
% subplot2d(nR,nC, [3 4], [3 4])

fakeTraceIndex = 4;

dph2 = loglog(lfpAmpNoise_sigmas, couplingSNR.plaEst(fakeTraceIndex,:)', ...
             'linewidth', 1, 'color','k');
hold all
dph1 = loglog(lfpAmpNoise_sigmas, couplingSNR.gplaEst(fakeTraceIndex,:)', ...
             'linewidth', 3, 'color','k');

nEst = size(couplingSNR.gplaEst, 1);
estColors = cool(nEst);

kEst_start = 3

for kEst = kEst_start : nEst
    loglog(lfpAmpNoise_sigmas, couplingSNR.plaEst(kEst, :), ...
                'linewidth', 1, 'color', estColors(kEst, :));
    hold on
end 
ax = gca;
ax.ColorOrderIndex = 1;

clear ph

for kEst = kEst_start : nEst
    ph(kEst-kEst_start+1) = loglog(lfpAmpNoise_sigmas, couplingSNR.gplaEst(kEst, :), ...
                'linewidth', 3, 'color', estColors(kEst, :));
    hold on
end
xlabel('LFP noise standard deviation \sigma')
% ylabel('Coupling SNR')
grid on

FRfacttorsTrimmed = FRfacttors(kEst_start:end);

% legend(ph, strcat('Rate factor=', cellstr(num2str(FRfacttors'))),'location','northeast')
[lh1, icons1] = legend(ph, strcat('Scale factor = ', ...
                                  cellstr(num2str(FRfacttorsTrimmed'))), 'location','northeast') 
% set(lh1, 'fontsize', vc.f9.gfs)
set(lh1, 'fontsize', vc.f9.lfs)

%%
bof = .4
% displaceFigureStuff(lh1, [.045 .04 NaN NaN])
% displaceFigureStuff(lh1, [.065 .04 NaN NaN])

% displaceFigureStuff(lh1, [.065 .04 0 0])
% old vals
% displaceFigureStuff(lh1, [.065 .05 0 0])
% new vals
displaceFigureStuff(lh1, [.065 .045 0 0])

for k = 6 : 2 : 15
    icons1(k).XData = icons1(k).XData + .3;
    tmpOrigPos = icons1(k).XData(2);
    % icons1(k).XData = tmpOrigPos + bof; 
    icons1(k).XData(1) = .05;
    icons1(k).XData(2) = icons1(k).XData(1) + .1
    % icons1(k).XData(2) = tmpOrigPos - bof; 
    icons1(k).LineWidth = 2;
end

for k = 1 : 7 - kEst_start+1
    % icons1(k).Position(1) = 0.4;
    icons1(k).Position(1) = 0.2;
    % icons1(k).FontSize = vc.f9.gfs
end

%%
% displaceFigureStuff(lh1, [.045 .04 NaN NaN])
% set(lh1,'FontSize', vc.f9.couplingCmpr_fs);
% set(gca, 'fontsize', vc.f9.gfs)
set(gca, 'fontsize', vc.f9.lfs)
set(lh1,'box', 'on');
lth1 = legendTitle(lh1, 'Color');
lth1.FontWeight = 'normal'
lth1.Position(2) = 0
lth1.Position(1) = .885;

% title(lh1, 'Color');
% lh1.Title.NodeChildren.Position = [0.1 .1 0];


% set(gca, 'fontsize', vc.f9.couplingCmpr_fs)
% set(gca, 'fontsize', vc.f9.gfs)
set_ticksOutward1(vc.f9.couplingCmpr_ts)
box off

ah1 = axes('position',get(gca,'position'),'visible','off');

[lh2, icons2] = legend(ah1, [dph1 dph2], 'GPLA-based', 'PLA-based', 'NorthWest');

% displaceFigureStuff(lh2, [NaN -.21 NaN NaN])
relocateFigureStuff(lh2, ...
                    [lh1.Position(1) + lh1.Position(3) - lh2.Position(3) ...
                    NaN NaN NaN])
% displaceFigureStuff(lh2, [NaN -.17 NaN NaN])
% old vals
% displaceFigureStuff(lh2, [NaN -.16 NaN NaN])

% new vals
displaceFigureStuff(lh2, [NaN -.12 NaN NaN])


lth2 = legendTitle(lh2, 'Thickness');
% set(lh2,'FontSize', vc.f9.couplingCmpr_fs);
set(lh2,'FontSize', vc.f9.gfs);
% set(lh2,'FontSize', vc.f9.lfs);

lth2.FontWeight = 'normal'
lth2.Position(2) = 0
lth2.Position(1) = .755;
% lth2.Position(1) = 0;