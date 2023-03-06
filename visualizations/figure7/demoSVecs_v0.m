selLFP = 3;
klfp = selLFP;


for kModel = 1 : 2
    switch kModel
      case 1
        coi = 0;  % col onset index - define which side of the figure it should locate

        load(fullfile(pds.datsrv.viz, 'neuron2020_submission',['neuralWaveStats1JLargeCenter' lfpExt{klfp}]), ...
             'lfpVec','spkVec','gPLV','allStat','map','globalDynamicsParams','SF','Nbin')

      case 2
        coi = 3;  % col onset index - define which side of the figure it should locate
        
        load(fullfile(pds.datsrv.viz, 'neuron2020_submission',['neuralWaveStats1SFJLargeCenter' lfpExt{klfp}]), ...
             'lfpVec','spkVec','gPLV','allStat','map','globalDynamicsParams','SF','Nbin')
        
    end


    % coi = 0;  % col onset index - define which side of the figure it should locate

    %%
    iofs = 4; % image offset
    % for kfreq = 1 : nfreq
    for kfreq = nfreq : -1 : 1

        % subplot2d(nR,nC, 4, (1:nhc))

        % remove zero values that bias the phase and replace by isotropic noise
        idxZeros = lfpVec.waveIllu{kfreq}==0;
        lfpVec.waveIllu{kfreq}(idxZeros) =  .0001*max(abs(lfpVec.waveIllu{kfreq}))*(randn(sum(idxZeros),1)+1i*randn(sum(idxZeros),1))';

        % cntr = cntr + 1;
        
        
        tri = 1 + coi; % temporary row index
        [~, sbh] = subplot2d(nR,nC, iofs + kfreq, tri);
        displaceFigureStuff(sbh, [hos(kModel) vos_svs NaN NaN])
        image_mapArray((spkVec.waveIllu{kfreq}), map);
        colormap hsv
        caxis([-pi pi])
        axis on; box on
        if kfreq == 1, 
            th = title('Spike vector'); 
            % th = title('Spike'); 
            set(th, 'fontweight', 'normal')
        end
        if coi == 0
            ylabel(sprintf('%d-%d Hz', globalDynamicsParams(kfreq).oscFreq))
        end
        set(gca,'xtick',[],'ytick',[])
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
        
        tri = 2 + coi; % temporary row index
        [~, sbh] = subplot2d(nR,nC, iofs + kfreq, tri);
        displaceFigureStuff(sbh, [hos(kModel) vos_svs NaN NaN])
        image_mapArray((lfpVec.waveIllu{kfreq}), map);
        colormap hsv
        caxis([-pi pi])
        axis on; box on
        % if cntr == 2, title('LFP vector'); end
        set(gca,'xtick',[],'ytick',[])
        if kfreq == 1, 
            th = title('LFP vector'); 
            % th = title('LFP'); 
            set(th, 'fontweight', 'normal')
            % xlabel('LFP vector')
        end
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0))
        
        
        % cntr = cntr + 1;
        % subplot(nR,nC, cntr);
        tri = 3 + coi; % temporary row index
        [~, sbh] = subplot2d(nR,nC, iofs + kfreq, tri);
        displaceFigureStuff(sbh, [hos(kModel) vos_svs NaN NaN])
        image_mapArray((spkVec.waveIllu{kfreq}.*conj(lfpVec.waveIllu{kfreq})), map);
        colormap hsv
        caxis([-pi pi])
        axis on; box on
        % if cntr == 3, title('Phase difference'); end
        set(gca,'xtick',[],'ytick',[])
        if kfreq == 1, 
            % th = title('Phase difference'); 
            th = title('Phase diff.'); 
            set(th, 'fontweight', 'normal')
        end
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0))

    end
end