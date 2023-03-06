%% LFP vector

% clf

freqBandsLab = {...
    '1-5';...
    '5-20'; ...
    '20-40'; ...
    '40-80'; ...
    '80-180' ...
    };

teh = 0.08; % temp enlarge height

[~, sbh] = subplot2d(nR,nC, [sdi+1 sdi+2], 1);
% displaceFigureStuff(sbh, [NaN NaN -.06 teh])
displaceFigureStuff(sbh, [NaN -.01 -.06 teh])

%%
% clf
maxElec = 32;
frCntr = 0;
for iFreq = foi
    frCntr = frCntr + 1;
    tsv = lfpVec.(caseName)(:,iFreq, isv);
    tlfpVecs(:,frCntr) = real(tsv);
    plot(real(tsv), (1:maxElec), 'linewidth',2)
    hold on
end

tmeanLfpVec = mean(tlfpVecs, 2);

[Electr,Ve] = compelectrode_two(ones(100,1),-ones(100,1),ones(100,1));
% the negative sign of Ve is just a scaling factor
% plot(-Ve, (1:maxElec), 'k:', 'linewidth',2) 

% we scale according to the maximum value
scaled_Ve = max(tmeanLfpVec) * -Ve / max(-Ve);
plot(scaled_Ve, (1:maxElec), 'k:', 'linewidth',2)

%%
% yline(8)
% yline(7)

ylabel('Channel ID')
xlabel('Magnitude')
set(gca, 'fontSize', vc.f11.gPLVxTickFS)

% [lh, icons] = legend(freqBandsLab(foi), 'location','southwest');
% set(lh, 'box', 'off')
% displaceFigureStuff(lh, [NaN .045 0 0])

allLabels = freqBandsLab(foi);
allLabels{numel(foi)+1} = 'Ground truth';
[lh, icons] = legend(allLabels, 'location','northwest');
set(lh, 'box', 'off')
% displaceFigureStuff(lh, [.03 -.01 0 0])
displaceFigureStuff(lh, [-.1 .012 0 0])

axis tight
lvl = .5; % LFP vector lim
xlim(lvl * [-1 1])
% xlim([-.5 0])

% set(gca, 'TickDir', 'out');
set(gca, ...
    'TickLength', [0.03 0.035], ...
    'TickDir',    'out', ...
    'ytick',      (8:8:32), ...
    'xtick',      [-lvl 0 lvl] ...
    );

% set(gca, 'TickDir', 'out');
% set(gca, 'TickLength', [0.0155 0.045]);
% set(gca, 'TickLength', [0.03 0.035]);

% grid on

for k = 1 : 4
    icons(k).Position(1) = 0.07;
end

bof = 0.3; % backward offset
for k = 5 : 2: 12
    tmpOrigPos = icons(k).XData(2);
    icons(k).XData(2) = tmpOrigPos - bof; 
end

for k = 5 : 2: 12
    tmpOrigPos = icons(k).XData(2);
    icons(k).XData(2) = tmpOrigPos - bof; 
end

ax = gca;
ax.YDir = 'reverse'
% ax.XDir = 'reverse'
box off

% % multi comparmtent model
% get the pos for multicomparmtn stuff, it willl plot later
posLfpVec = sbh.Position;
% axes('Position',[posLfpVec(1)+posLfpVec(3)+.01 ...
%                  posLfpVec(2)-.001 .1 posLfpVec(4)*1.05])
% imshow('multCompModelNeuron_v1_crop.png')
% set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

% % neuron multi compartment
% teh = 0.02; % temp enlarge height
% clf

tmos = .005; % additional offset to bring the soma of the
             % multicomparmet model with elec 7

[~, sbh] = subplot2d(nR,nC, [sdi+1 sdi+2], 2);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
displaceFigureStuff(sbh, [NaN -tmos NaN NaN]);
% teh = 0.08; % temp enlarge height
displaceFigureStuff(sbh, [-.049 -.005 .1 teh]);
imshow('detailedModel_v2.png')
% axis on
% %
% placeHolder
% subplot2d(nR,nC, [sdi+2 sdi+3], 3)
% box on

% % multi comparmtent model - position was defined earlier (close
% to LFP )
% axes('Position',[posLfpVec(1)+posLfpVec(3)+.001 ...
%                  posLfpVec(2)-.001 .1 posLfpVec(4)*1.05])
% axes('Position',[posLfpVec(1)+posLfpVec(3)-.02 ...
%                  posLfpVec(2)-.001 .1 posLfpVec(4)*1.05])

axes('Position',[posLfpVec(1)+posLfpVec(3)-.02 ...
                 posLfpVec(2)-.001-tmos .1 posLfpVec(4)*1.05])
imshow('multCompModelNeuron_v1_crop.png')
set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

%% compute cosine similarity and store
frCntr = 0;
clear cosSim;
for iFreq = foi

    frCntr = frCntr + 1;
    tsv = real(lfpVec.(caseName)(:,iFreq, isv));
    cosSimInfo(frCntr).sim = cmpt_cosineSimilarity(tsv, scaled_Ve);
    cosSimInfo(frCntr).freqLab = freqBandsLab{iFreq};
    
end

cosSimInfoTable = struct2table(cosSimInfo);
% *** save it

%% draft - ground truth
% clf
% maxElec = 32;
% for iFreq = foi
%     frCntr = frCntr + 1;
%     tsv = lfpVec.(caseName)(:,iFreq, isv);
%     plot(real(tsv), (1:maxElec), 'linewidth',2)
%     hold on
% end


% % [Electr,Ve] = compelectrode_two(ones(100,1),-ones(100,1),ones(100,1));
% % plot(-Ve, 2*Electr(:,:,2) / maxElec, 'k--') 
% plot(-Ve, (1:maxElec), 'k--') 
% % set(gca,'ydir','reverse')

