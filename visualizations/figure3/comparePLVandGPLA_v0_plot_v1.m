% this slightly different version of comparePLVandGPLA_v0_plot

%%
barData(1 : signalParams.nUnit) = 100 * sum(statSummary.PLV) / nRel;
barData(signalParams.nUnit + 1) = 100 * sum(statSummary.pPLV) / nRel;
barData(signalParams.nUnit + 2) = 100 * sum(statSummary.gPLV) / nRel;


%%
[~,~,~,~,~, rawSvdStuff.(caseType)(iCase)] = ...
        tngpla(spikeTrains.(caseType){iCase}, lfpLikePhase.(caseType){iCase}, [], [], [], [] , [], ...
               [], iSV, sameElecCheckInfo, 'nSpk', 0);

subplot2(nR,nC, 2, 1);

for iRel = 1 : nRel
    tmp = PLVexamples(iRel, :);
    polarplot(angle(tmp), abs(tmp),...
              '.', 'color',.6*ones(1,3), 'MarkerSize', vc.f9.bulletSize_tiny); 
    hold on
end


% tmp = rawSvdStuff.gplaBenefit.couplingMatrix;
%%
vc.f9.bulletSize = 10;
tmp = mean(PLVexamples, 1);
polarplot(angle(tmp), abs(tmp),...
    'k.', 'MarkerSize', vc.f9.bulletSize);
%%
thetaticks(0:45:315)
tn = .5;
rlim([0 tn])
rticks([0 tn/2 tn]);


ax = gca;
ax.LineWidth = 1.5;
rruler = ax.RAxis;
rruler.Label.String = 'PLV';
ax.RAxisLocation = 45/2;
set(gca, 'fontsize', vc.f9.gfs)


%%

subplot2(nR,nC, 2, 2);
bar([1 2 3 5 7], barData, vc.f7.barWidth, 'facecolor', 'k')
grid on

ylabel('Detection ratio [%]')

% xticklabels({'Unit 1 PLV','Unit 2 PLV','Unit 3 PLV', 'pPLV', 'gPLV'}); %
% xticklabels({'Unit 1','Unit 2','Unit 3', 'pPLV', 'gPLV'});
xticklabels({'Unit 1','Unit 2','Unit 3', 'pPLV', 'gPLV'});
% xticklabels({'','PLV','Unit 3', 'pPLV', 'gPLV'});
xtickangle(45)

box off
mr = 100; % max ratio
% ylim([0 mr])
ylim([0 mr])
set(gca, 'ytick', (0 : 20 : 100))

set_ticksOutward1(.05)

% will moved out
vc.f9.gfs = 6; % global font size
set(gca, 'fontsize', vc.f9.gfs)
