%%
ncc = 1;  % we only play with one case of cell clusterting
iCase = 1; % is useless here

allOscFreq = 11:.5: 15.5;
globalDynamicsParams.oscFreq = allOscFreq;
globalDynamicsParams.nFreqComp = numel(globalDynamicsParams.oscFreq);

%%
noc = numel(globalDynamicsParams.oscFreq);
% *** should be removed
% oscWeights are useless, given we have mixing matrix here, for
% sake syntax compatibilty  we keep them
globalDynamicsParams.oscComps = ones(1, noc);


nclu = globalDynamicsParams.nFreqComp;

%%
iun = 1;
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

%%
% ~ allocate 
tmpSummaryStat = NaN(nPopNum, nCouplingStrength, nRel);

for ipn = 1 : nPopNum
    
    cer = eye(nclu); % clusterer: define the population structure
    cer((ipn+1:nclu), (ipn+1:nclu)) = 0;

    % mixing matrix
    baseLineCoef = 0.1; % 0 correspond to no mixing
    globalDynamicsParams.mixingMatrix = kron(eye(noc), (1 - baseLineCoef)*ones(signalParams.nCh/noc,1)) + baseLineCoef;

    baseLineCoef = 0;
    spikeTrainParams.avefiringRate = kron(eye(noc), aveFR * (1 - baseLineCoef)*ones(signalParams.nUnit/noc,1)) + baseLineCoef;

    for ic = 1 : nCouplingStrength
        ic
        icc = 1; % usless index, should be removed

        couplingStrength = couplingStrengths(ic);
        % coupling strength is specified for all population
        baseLineCoef = 0; % 0 correspond to no coupling across other freq
        couplingParams(iCase).lockingStrength_kappa = ...
            kron(cer(:,:, icc), couplingStrength * (1 - baseLineCoef)*ones(signalParams.nUnit/ noc, 1)) + baseLineCoef;
        
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
            [~,~,~,~, allStats, svdOut] = ...
                tngpla(spikeTrains, analLfp, [], [], [], [] , [], ...
                       statTestInfo, 'all', sameElecCheckInfo, ...
                       'nSpk-square-root', 2, 0);

            % tmpSummaryStat(ipn, ic, iRel) = allStats.gPLV_stats.nullHypoReject;
            tmpSummaryStat(ipn, ic, iRel) = sum(allStats.SV_stats.nullHypoReject);
        end
    end
end                                  
summaryStat.(caseName) = tmpSummaryStat;


%% tests
% tmpSummaryStat = NaN(nPopNum, nCouplingStrength, nRel);

% ipn = 7

% cer = eye(nclu); % clusterer: define the population structure
% cer((ipn+1:nclu), (ipn+1:nclu)) = 0;

% % mixing matrix
% baseLineCoef = 0.1; % 0 correspond to no mixing
% globalDynamicsParams.mixingMatrix = kron(eye(noc), (1 - baseLineCoef)*ones(signalParams.nCh/noc,1)) + baseLineCoef;

% baseLineCoef = 0;
% spikeTrainParams.avefiringRate = kron(eye(noc), aveFR * (1 - baseLineCoef)*ones(signalParams.nUnit/noc,1)) + baseLineCoef;

% ic = 10
% icc = 1; % usless index, should be removed

% couplingStrength = couplingStrengths(ic);
% % coupling strength is specified for all population
% baseLineCoef = 0; % 0 correspond to no coupling across other freq
% couplingParams(iCase).lockingStrength_kappa = ...
%     kron(cer(:,:, icc), couplingStrength * (1 - baseLineCoef)*ones(signalParams.nUnit/ noc, 1)) + baseLineCoef;


% parfor iRel = 1 : nRel
%     %% simulation 
%     [lfpLikeSig, spikeTrains] = ...
%         smlt_sustLockedSpkLfpPairs_multFreqCoupling2(...
%             globalDynamicsParams, spikeTrainParams, couplingParams, ...
%             signalParams);

%     %% filter LFP 
%     [~, analLfp] = tpp_filt_recenter(lfpLikeSig, freqCenter, halfFilterWidth, ...
%                                      signalParams.SF, filterOrder, nIteCent);

%     %% GPLA
%     [~,~,~,~, allStats, svdOut] = ...
%         tngpla(spikeTrains, analLfp, [], [], [], [] , [], ...
%                statTestInfo, 'all', sameElecCheckInfo, ...
%                'nSpk-square-root', 2, 0);

%     % tmpSummaryStat(ipn, ic, iRel) = allStats.gPLV_stats.nullHypoReject;
%     tmpSummaryStat(ipn, ic, iRel) = sum(allStats.SV_stats.nullHypoReject);
% end

% ipn
% ic
% squeeze(tmpSummaryStat(ipn,ic,:))