%%
ncc = 1;  % we only play with one case of cell clusterting

nclu = globalDynamicsParams.nFreqComp;

%%
iCase = 1; % is useless here

% ~ allocate 
tmpSummaryStat = NaN(nUnitNum, nRel);

for iun = 1 : nUnitNum
    signalParams.nCh = unitNums(iun);
    signalParams.nUnit = unitNums(iun);

    % coupling strength is specified later for the population
    couplingParams(iCase) = struct ...
        (...
            'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
            0, ...
            'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
            2*pi * rand(signalParams.nUnit, 1) ... % rad  
            );

    
    % mixing matrix
    baseLineCoef = 0.1; % 0 correspond to no mixing
    globalDynamicsParams.mixingMatrix = kron(eye(noc), (1 - baseLineCoef)*ones(signalParams.nCh/noc,1)) + baseLineCoef;

    % coupling strength is specified for all population
    couplingParams(iCase).lockingStrength_kappa = zeros(signalParams.nUnit, globalDynamicsParams.nFreqComp);

    baseLineCoef = 0;
    spikeTrainParams.avefiringRate = kron(eye(noc), aveFR * (1 - baseLineCoef)*ones(signalParams.nUnit/noc,1)) + baseLineCoef;

    parfor iRel = 1 : nRel
        iRel
        %% simulation 
        [lfpLikeSig, spikeTrains] = ...
            smlt_sustLockedSpkLfpPairs_multFreqCoupling2(...
                globalDynamicsParams, spikeTrainParams, couplingParams, ...
                signalParams);

        %% filter LFP 
        [~, analLfp] = tpp_filt_recenter(lfpLikeSig, freqCenter, halfFilterWidth, ...
                                                    signalParams.SF, filterOrder, nIteCent);

        %% GPLA
        icc = 1; % usless index, should be removed

        [~,~,~,~, allStats, svdOut] = ...
            tngpla(spikeTrains, analLfp, [], [], [], [] , [], ...
                   statTestInfo, 'all', sameElecCheckInfo, ...
                   'nSpk-square-root', 2, 0);

        tmpSummaryStat(iun, iRel) = allStats.gPLV_stats.nullHypoReject;
    end
end                                  
summaryStat.(caseName) = tmpSummaryStat;