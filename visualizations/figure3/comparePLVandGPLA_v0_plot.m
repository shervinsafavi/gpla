%%
barData(1 : signalParams.nUnit) = 100 * sum(statSummary.PLV) / nRel;
barData(signalParams.nUnit + 1) = 100 * sum(statSummary.pPLV) / nRel;
barData(signalParams.nUnit + 2) = 100 * sum(statSummary.gPLV) / nRel;


%%
[~,~,~,~,~, rawSvdStuff.(caseType)(iCase)] = ...
        tngpla(spikeTrains.(caseType){iCase}, lfpLikePhase.(caseType){iCase}, [], [], [], [] , [], ...
               [], iSV, sameElecCheckInfo, 'nSpk', 0);

subplot2(nR,nC, 1, 3);

for iRel = 1 : nRel
    tmp = PLVexamples(iRel, :);
    polarplot(angle(tmp), abs(tmp),...
              '.', 'color',.6*ones(1,3), 'MarkerSize', vc.f7.bulletSize/5); 
    hold on
end


% tmp = rawSvdStuff.gplaBenefit.couplingMatrix;
tmp = mean(PLVexamples, 1);
polarplot(angle(tmp), abs(tmp),...
    'k.', 'MarkerSize', vc.f7.bulletSize);

rlim([0 .5])
ax = gca;
ax.LineWidth = 1.5;
rruler = ax.RAxis;
rruler.Label.String = 'PLV';
ax.RAxisLocation = 45;

%%

subplot2(nR,nC, 1, 4);
bar(barData, vc.f7.barWidth, 'facecolor', 'k')
grid on

ylabel('Detection ratio [%]')

% xticklabels({'Unit 1 PLV','Unit 2 PLV','Unit 3 PLV', 'pPLV', 'gPLV'}); %
% xticklabels({'Unit 1','Unit 2','Unit 3', 'pPLV', 'gPLV'});
xticklabels({'Unit 1','Unit 2','Unit 3', 'pPLV', 'gPLV'});
% xticklabels({'','PLV','Unit 3', 'pPLV', 'gPLV'});
xtickangle(45)

box off
mr = 60; % max ratio
ylim([0 mr])
set(gca, 'ytick', [0 mr/2 mr])

set_ticksOutward1(.05)

% will moved out
% vc.f9.gfs = 6; % global font size
set(gca, 'fontsize', vc.f9.gfs)
