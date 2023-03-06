% cer = []; % clusterer: define the population structure
ncc = 1;  % we only play with one case of cell clusterting
% nclu = globalDynamicsParams.nFreqComp; % number of clusters (because we deal with 6 freq)
                                       % cer = eye(nclu);
                                       % cer((3:6), (3:6)) = 0;

nclu = globalDynamicsParams.nFreqComp;

%
iCase = 1; % is useless here

% ~ allocate 
% summaryStat.(caseName) = NaN(nUnitNum, nCouplingStrength, nRel);
tmpSummaryStat = NaN(nUnitNum, nCouplingStrength, nRel);

for iun = 1 : nUnitNum
    iun
    signalParams.nCh = unitNums(iun);
    signalParams.nUnit = unitNums(iun);

    % coupling strength is specified later for the population
    couplingParams(iCase) = struct ...
        (...
            'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
            NaN, ...
            'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
            2*pi * rand(signalParams.nUnit, 1) ... % rad  
            );

    % for icc = 1 : ncc % loop is not necessary here
    % globalDynamicsParams.oscFreq = allOscFreq;
    
    % mixing matrix
    baseLineCoef = 0.1; % 0 correspond to no mixing
    globalDynamicsParams.mixingMatrix = kron(eye(noc), (1 - baseLineCoef)*ones(signalParams.nCh/noc,1)) + baseLineCoef;

    baseLineCoef = 0;
    spikeTrainParams.avefiringRate = kron(eye(noc), aveFR * (1 - baseLineCoef)*ones(signalParams.nUnit/noc,1)) + baseLineCoef;

    for ic = 1 : nCouplingStrength
        ic
        couplingStrength = couplingStrengths(ic);
        % coupling strength is specified for all population
        couplingParams(iCase).lockingStrength_kappa = zeros(signalParams.nUnit, globalDynamicsParams.nFreqComp);
        couplingParams(iCase).lockingStrength_kappa(:, 1) = couplingStrength * ones(signalParams.nUnit, 1);


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
            tmpSummaryStat(iun, ic, iRel) = allStats.gPLV_stats.nullHypoReject;
        end
    end
end                                     %
summaryStat.(caseName) = tmpSummaryStat;
