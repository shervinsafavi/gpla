% kModel = 1;

for kModel = 1 : 2
    switch kModel
      case 1
        coi = 0;  % col onset index - define which side of the figure it should locate

        % *** change to the currect path
        % fileName = 'simneuralwaveJLarge_28265.mat';
        % fileNameSpk = 'spkHR_simneuralwaveJLarge_28265.mat';
        fileName = fullfile(pds.datsrv.viz, 'neuron2020_submission', 'simneuralwaveJLarge_28265.mat');
        fileNameSpk = fullfile(pds.datsrv.viz, 'neuron2020_submission', 'spkHR_simneuralwaveJLarge_28265.mat');
      case 2
        coi = 3;  % col onset index - define which side of the figure it should locate
        
        fileName = fullfile(pds.datsrv.viz, 'neuron2020_submission', 'simneuralwaveSFJLarge_13062.mat');
        fileNameSpk = fullfile(pds.datsrv.viz, 'neuron2020_submission', 'spkHR_simneuralwaveSFJLarge_13062.mat');
    end
    
    if ~exist('idxRange','var')
        idxRange = 1:length(lfpIdx);
    end

    klfp = idxRange;

    clear spikeTrainTot lfpLikeAnalSig
    iCase = 1;%:length(globalDynamicsParams)
    freqBand = globalDynamicsParams(iCase).oscFreq;
    
    %% statistics parameters 
    clear statTestInfo
    statTestInfo.testType = 'spike-jittering';
    statTestInfo.SVspectrumStatsType = 'RMT-heuristic';
    statTestInfo.jitterType = 'interval-jittering';

    statTestInfo.nJtr = 100;
    statTestInfo.alphaValue = 0.05;
    
    statTestInfo.jitterWinWidth = 1 / mean(freqBand);




    clear lfpLikeAnalSig

    iTr = 1;%1 : min(50,length(fileList))
            %res = load([fileList(iTr).folder filesep fileList(iTr).name],'SimRes','EdpParams','StimVec','SimVarLabels');
    res = load(fileName,'SimRes','EdpParams','StimVec','SimVarLabels');
    %spkTrain = load([fileList(iTr).folder filesep 'spkHR_' fileList(iTr).name],'spikeTrains');
    spkTrain = load(fileNameSpk,'spikeTrains');

    %EmeanRate = 5; % use 5Hz as typical pyramidal cell rate
    % filter

    SF = 1/(res.EdpParams.StoreSubSampT*res.EdpParams.dt);
    statTestInfo.spkSF = SF;

    filterOrder = 2;

    wn = freqBand / (SF/2);
    [b, a] = butter(filterOrder, wn, 'bandpass');
    
    
    gridDims = [res.EdpParams.szx,res.EdpParams.szy]/res.EdpParams.StoreSubSampX/res.EdpParams.dx;
    spatialFreq = 1/res.EdpParams.StoreSubSampX/res.EdpParams.dx;

    Nbin = size(res.SimRes,3);
    %ErateResh = reshape(res.SimRes(:,:,:,5),[],Nbin);
    spikeTrainDuration  =  Nbin/SF;
    %spikeTrainE = gnrt_inhomogeneousPoissonSpikeTrains(ErateResh*EmeanRate, spikeTrainDuration, SF);

    centerIdx = 6:20;
    
    tidx = 1000:2000;
    % figure

    % [~, sbh] = subplot2d(nR,nC, 1, coi + nhc - 1)
    % displaceFigureStuff(sbh, [hos(kModel) NaN NaN NaN])    

    % load data script

    % input
    [~, sbh] = subplot2d(nR,nC, 2, coi + (1:nhc))
    displaceFigureStuff(sbh, [hos(kModel) NaN NaN NaN])
    
    % plot(squeeze(res.SimRes(15,15,tidx,[ end])),'color',[0,.5,0])
    plot(squeeze(res.SimRes(15,15,tidx,[ end])),'color','k')
    %legend(res.SimVarLabels([end]))
    % legend('Exogenous input')
    % tyl = sprintf('Post-synaptic \n current (a.u.)')
    % ylabel('Post-synaptic \n current (a.u.)')
    if kModel == 1
        ylabel({'Post-synaptic','current [a.u.]'})
    end
    nicefig
    axis tight

    %%
    [~, sbh] = subplot2d(nR,nC, 3, coi + (1:nhc))
    displaceFigureStuff(sbh, [hos(kModel) NaN NaN NaN])
    plot(20*squeeze(res.SimRes(15,15,tidx,5:end-1)))
    % legend(res.SimVarLabels(5:end-1))
    if kModel == 1
        lh = legend('Excitatory', 'Inhibitory')
        displaceFigureStuff(lh, [.07 -.008 nan nan])
    end
    nicefig
    axis tight
    if kModel == 1
        ylabel({'Population spike', 'rate [Hz]'})
    end
    %
    % clf
    
    %%
    map = reshape(1:size(spkTrain.spikeTrains,1),25,25);
    % spkTS = full(spkTrain.spikeTrains(map(15,15),:));
    allSpikes = full(spkTrain.spikeTrains(:,tidx));
    
    [~, sbh] = subplot2d(nR,nC, 4, coi + (1:nhc));
    % displaceFigureStuff(sbh, [hos(kModel) NaN NaN NaN])
    displaceFigureStuff(sbh, [hos(kModel) -.02 NaN NaN])
    % plot_eventRaster(allSpikes, {'color','k'});
    plot_eventRaster(allSpikes, {'color','k', 'markersize',2.5}, [], 'scatter');
    plotScaleBar1([0 size(allSpikes, 2)], [1 size(allSpikes, 1)], '200 ms', NaN, 2 * SF/10, NaN, ...
              -1.1, NaN, -.125, NaN, 1, 6.5); 
    box off
    axis off
    
end
% %% viz (only for one)
% % E-I schematics
% subplot2d(nR,nC, 1, coi + nhc - 1)

% % load data script

% % input
% subplot2d(nR,nC, 2, coi + (1:nhc))
% % plot(squeeze(res.SimRes(15,15,tidx,[ end])),'color',[0,.5,0])
% plot(squeeze(res.SimRes(15,15,tidx,[ end])),'color','k')
% %legend(res.SimVarLabels([end]))
% % legend('Exogenous input')
% % tyl = sprintf('Post-synaptic \n current (a.u.)')
% % ylabel('Post-synaptic \n current (a.u.)')
% ylabel({'Post-synaptic','current [a.u.]'})

% nicefig
% axis tight


% subplot2d(nR,nC, 3, coi + (1:nhc))
% plot(20*squeeze(res.SimRes(15,15,tidx,5:end-1)))
% legend(res.SimVarLabels(5:end-1))
% nicefig
% axis tight
% ylabel({'Population spike', 'rate (Hz)'})

% %
% % clf
% subplot2d(nR,nC, 4, coi + (1:nhc));
% % plot_eventRaster(allSpikes, {'color','k'});
% plot_eventRaster(allSpikes, {'color','k', 'markersize',2.5}, [], 'scatter');
% box off
% axis off
