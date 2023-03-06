% ignit_gpla
% ignit_datsrv
% pds = ignit_datsrv(pds)

%%
lfpIdx = [3,7,1,2];
lfpExt = {'Ve','exo','epsp','ipsp'};
simList = [1:2];
simLabel = {'Weak recurrence','Strong recurrence'} ;
selBand =1:4;%1:length(globalDynamicsParams)
selLFP = 3;

plotTypeIds = [1 4 3];

% figure
for klfp = selLFP
    clear maxPLVfreq
    kplotidx = 0;
    for kPlt = 1 : numel(plotTypeIds)
        kplot = plotTypeIds(kPlt);
        % kplotidx =kplotidx+1;
        % subplot(1,3,kplotidx)
        [~, sbh] = subplot2d(nR, nC, trInd, tcInd(kPlt,:));
        if kPlt == 1
            displaceFigureStuff(sbh, [-hos2 vos_svs2 NaN NaN])
        elseif kPlt == 2
            displaceFigureStuff(sbh, [NaN vos_svs2 NaN NaN])
        elseif kPlt == 3
            displaceFigureStuff(sbh, [+hos2 vos_svs2 NaN NaN])
        end
        for ksim = simList
            % subplot2d(2,4, ksim,kplotidx+1)
            

            switch ksim
                
              case 1
                
                load(fullfile(pds.datsrv.viz, 'neuron2020_submission',['neuralWaveStats1JLargeCenter' lfpExt{klfp}]),'lfpVec','spkVec','gPLV','allStat','map','globalDynamicsParams','SF','Nbin')

              case 2
                load(fullfile(pds.datsrv.viz, 'neuron2020_submission', ...
                              ['neuralWaveStats1SFJLargeCenter' lfpExt{klfp}]),'lfpVec','spkVec','gPLV','allStat','map','globalDynamicsParams','SF','Nbin')
            end
            
            [distances, distancesTabel] = cmpt_arrayPWdistancesFromCenter(map);
            bands = cat(1,globalDynamicsParams(:).oscFreq);

            switch kplot
              case 1
                
                plvVals =    [gPLV.waveIllu{:}];
                [~,maxPLVfreq(ksim)] = max(plvVals);
                plot(mean(bands,2),plvVals, 'linewidth', 2, 'color', ...
                     netColors(ksim, :))

                box off
                set(gca, 'TickDir', 'out');
                set(gca, 'TickLength', [0.0355 0.035]);

                lh = legend('Weak recurrence', 'Strong recurrence')
                set(lh, 'box', 'off')
                displaceFigureStuff(lh, [.02 .001 NaN NaN])
                % displaceFigureStuff(lh, [.12 .0055 NaN NaN])
                
                hold on

                % legend(simLabel(simList(1:ksim)))
                
                xlabel('Frequency [Hz]')
                ylabel('gPLV')
                % title('gPLV')
              case 2
                % restrict plot to 6mm radius around center
                semilogx(abs(lfpVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)), 180/pi*angle(lfpVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)),'.');
	        
                hold on
                xlabel('LFP vector coef. modulus (a.u.)')
                ylabel('LFP vector coef. phase (deg.)')
                title('Spatial variations')
              case 3
                
                plot(mean(bands,2),180/pi*angle(cellfun(@mean, ...
                                                        spkVec ...
                                                        .waveIllu)),'Linewidth',2, 'color', netColors(ksim, :))
                box off
                set(gca, 'TickDir', 'out');
                set(gca, 'TickLength', [0.0355 0.035]);
                
                hold on
                xlabel('frequency [Hz]')
                ylabel('Phase diff. [deg]')
                % title('spike-LFP shift')
              case 4
                % restrict plot to 6mm radius around center
                h1 =plot(abs(spkVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)), 180/pi*angle(spkVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)),'.', 'color', netColors(ksim, :));
                %plot(abs(spkVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)), 180/pi*angle(lfpVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)),'.');
                res = fitlm((abs(spkVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4))),180/pi*angle(spkVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)));
                
                hold on
                x =(abs(spkVec.waveIllu{maxPLVfreq(ksim)}(distances*.8<4)));
               
                box off
                set(gca, 'TickDir', 'out');
                set(gca, 'TickLength', [0.0355 0.035]);                
               
                ylim([-50 100])
                plot(x,x*res.Coefficients.Estimate(2)+ ...
                     res.Coefficients.Estimate(1),'color', ...
                     h1.Color,'Linewidth',2, 'color', netColors(ksim, :))
                % xlabel('Spike vector coef. modulus (a.u.)')
                % ylabel('Spike vector coef. phase (deg.)')
                xlabel('Modulus [a.u.]')
                ylabel('Phase [deg]')
                % title('Spatial variations')
                
            end
            
        end
        % nicefig
    end
end

%%
% eiNetImages = {'EInet_weakRecurrent.png', ...
%                'EInet_storngRecurrent.png'};

% nModel = 2;

% % vos_eiNets = .03 - 0.02;
% vos_eiNets = .01;
% hos_eiNets = .04;

% for kModel = 1 : nModel
%     % if kModel == 1, coi = 0; 
%     % elseif kModel == 2, coi = 3; end
%     % [~, sbh] = subplot2d(nR,nC, 1, coi + nhc - 1)
%     % displaceFigureStuff(sbh, [hos(kModel)-hos_eiNets vos_eiNets .1 .01])    
%     % imshow('EInet_weakRecurrent.png');

%     subplot2d(2,4, kModel, 1);
    
%     [A, map, transparency] = imread(eiNetImages{kModel});
%     image(A, 'AlphaData' , transparency);
%     % imshow (A, map, 'Parent' , gc)
%     % imshow (A)
%     axis off
%     pbaspect([1 .65 1])
%     % set(gca, 'Color' , [1 1 1] )
% end

%%
% [~, sbh] = subplot2d(2,4, 1, 4);
% displaceFigureStuff(sbh, [.035 -.25 NaN NaN])    


%%
% % terebele naming
% fsp = '/home/ssafavi/Nextcloud/research/nnr/reports_nnr/papers/p_ploscb2018_gpla/figures'
% svflg = 1;
% fn = 'abm2020_fig1_neualFieldDemo';
% spf(fsp, fn, svflg)
