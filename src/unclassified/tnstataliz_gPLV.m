function [gPLV, pValue,  lfpVec, spkVec, couplingMatrix, singularValues, ...
          varargout] = tnstataliz_gPLV(spikeTrains, lfpPhases, ...
                                       statTestInfo, iSV, ...
                                       sameElecCheckInfo, plvNrmlzMethed, ...
                                       unwhitenOpr, flag_gPLVnrmlz)
% [gPLV, pValue,  lfpVec, spkVec,couplingMatrix, singularValues, nullHypoReject,
%       gPLV_stats, PLV_stats, SVspectrum_stats] = stataliz_gPLV(spikeTrains, lfpPhases, statTestInfo, iSV)
%
% gpla significance testbased on surrogate data for null hypothesis.
% EXAMPLE:
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     spikeTrains: 
%       expnation
% 2     lfpPhase: 
%       expnation
% 3   statTestInfo: structure
%       (default value = 1)
%       it varies according to method, for:
%       'permutation'
%           testType
%           nPrm
%           alphaValue
%       number of cells (i.e. generate Poisson spike rrain for a population
% 4     selected singular value (default = 1)
%
% Output:
% 1     gPLV, 
% 2     pValue,
% 3     lfpVec
% 4     spkVec
% 5     couplingMatrix
% 6     singularValues
% 
% (7)   nullHypoReject: logical outcome of test (True if rejection)
% (8)   gPLV_stats:
%          CI: 
%          nullDistribution: 
% (9)   PLV_stats: 
%          CI: 
%          nullDistribution: 
%          pValue: 
% (10) SVspectrum_stats
% ------
% potential improvments:
% (1) flag for online updates on the progess (e.g. iJtr and pVal)
% (2) handling the case where uder passed one dimentional stuff
% (specially for p-value calcaulation
% ------
% Code Info:
%   creation: 2019-06-04 by SS and MB (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ 2018-??-?? ???
% ------
% see also
% output: null distribution for PLV and gPLV

%% Data info
n.LfpCh     = size(lfpPhases, 1);        % number of LFP channels
n.Sample    = size(lfpPhases, 2);        % number of samples
n.SpkUnit   = size(spikeTrains, 1);   % number of spiking units


if isfield(statTestInfo, 'nJtr')
    statTestInfo.nSd = statTestInfo.nJtr;
elseif isfield(statTestInfo, 'nPrm')
    statTestInfo.nSd = statTestInfo.nPrm;    
end

% statTestInfo.jitterType = 'interval-jittering';

% *** if sameElecCheckInfo is passed with some missing fields, fill it


%%
[lfpVec, spkVec, gPLV, cgPLV, couplingMatrix, singularValues] = ...
    gpla_core(spikeTrains, lfpPhases, 0, iSV, sameElecCheckInfo, ...
              plvNrmlzMethed, unwhitenOpr);

% if flag_gPLVnrmlz == 1
%     [lfpVec, spkVec, gPLV_out, cgPLV, couplingMatrix, singularValues] = ...
%         gpla_core(spikeTrains, lfpPhases, flag_gPLVnrmlz, iSV, sameElecCheckInfo, ...
%                   plvNrmlzMethed, unwhitenOpr);
%     % we need the non-normalized one for statistics 
%     [lfpVec, spkVec, gPLV, cgPLV, couplingMatrix, singularValues] = ...
%         gpla_core(spikeTrains, lfpPhases, 0, iSV, sameElecCheckInfo, ...
%                   plvNrmlzMethed, unwhitenOpr);
% else
%     [lfpVec, spkVec, gPLV, cgPLV, couplingMatrix, singularValues] = ...
%         gpla_core(spikeTrains, lfpPhases, 0, iSV, sameElecCheckInfo, ...
%                   plvNrmlzMethed, unwhitenOpr);
% end


% singularValues = singularValues'; % in some cases (1D) this is
                                  % better to avoide errors

singularValues = singularValues(:);

% *** not sure if that a proper solution
% if iSV > 1
%     singularValues = singularValues';
% end

% if isempty(sameElecCheckInfo)
%     [lfpVec, spkVec, gPLV, cgPLV, couplingMatrix, singularValues] = ...
%         gpla_core(spikeTrains, lfpPhases, 0, iSV);
% else 
%     [lfpVec, spkVec, gPLV, cgPLV, couplingMatrix, singularValues] = ...
%         gpla_core(spikeTrains, lfpPhases, 0, iSV, sameElecCheckInfo);
% end

%% run the statistics if applicapble

if ~isempty(statTestInfo)
    % choose the method
    switch statTestInfo.testType
        
      case 'spike-jittering'
        % ~ meomry allocation
        jtrdspt_gPLVs             = nan(statTestInfo.nJtr, 1); % jtrdspt: Jittered Spike Train
        jtrdspt_couplingMatrix    = nan(n.LfpCh, n.SpkUnit, statTestInfo.nJtr);
        %         jtrd_spikeTrains{1,n.Tr} = nan(size(spikeTrains{n.Tr}));
        fprintf('spike jittering iterations\n')
        % for iJtr = 1 : statTestInfo.nJtr
        %             iJtr
        % *** only works, if binarySparseSpikeTrainpassed to the function            
        switch statTestInfo.jitterType
          case 'interval-jittering'
            for iJtr = 1 : statTestInfo.nJtr
                jtrd_spikeTrains = ...
                    jitter_binarySparseSpikeTrain(spikeTrains, statTestInfo.jitterWinWidth, statTestInfo.spkSF);
                fprintf('\n')
                [~,~, jtrdspt_gPLVs(iJtr),  ~, jtrdspt_couplingMatrix(:,:, iJtr) , jtrdspt_singularValues(:, iJtr)] = ...
                    gpla_core(jtrd_spikeTrains, lfpPhases, 0, iSV, ...
                              [], plvNrmlzMethed, unwhitenOpr);
                fprintf('.')            
            end

          case 'ISI-preserved-interval-jittering'
            for iJtr = 1 : statTestInfo.nJtr
                jtrd_spikeTrains = ...
                    jitter_binarySparseSpikeTrainISI(spikeTrains, statTestInfo.jitterWinWidth, statTestInfo.spkSF);
                fprintf('\n')
                [~,~, jtrdspt_gPLVs(iJtr),  ~, jtrdspt_couplingMatrix(:,:, iJtr) , jtrdspt_singularValues(:, iJtr)] = ...
                    gpla_core(jtrd_spikeTrains, lfpPhases, 0, iSV, [], plvNrmlzMethed, unwhitenOpr);
                fprintf('.')            
            end
          case 'group-preserved-interval-jittering'
            for iJtr = 1 : statTestInfo.nJtr
                jtrd_spikeTrains = ...
                    jitter_binarySparseSpikeTrainGroup(spikeTrains, statTestInfo.jitterWinWidth, statTestInfo.spkSF);
                fprintf('\n')
                [~,~, jtrdspt_gPLVs(iJtr),  ~, jtrdspt_couplingMatrix(:,:, iJtr) , jtrdspt_singularValues(:, iJtr)] = ...
                    gpla_core(jtrd_spikeTrains, lfpPhases, 0, iSV, [], plvNrmlzMethed, unwhitenOpr);
                fprintf('.')            
            end
          case 'sync-jittering'
            for iJtr = 1 : statTestInfo.nJtr
                jtrd_spikeTrains = ...
                    popJitter_binarySparseSpikeTrain_multWW(spikeTrains, ...
                                                            statTestInfo.jitterWinWidthRange, statTestInfo.spkSF);
                [~,~, jtrdspt_gPLVs(iJtr),  ~, jtrdspt_couplingMatrix(:,:, iJtr) , jtrdspt_singularValues(:, iJtr)] = ...
                    gpla_core(jtrd_spikeTrains, lfpPhases, 0, iSV, [], plvNrmlzMethed, unwhitenOpr);
                fprintf('.')
            end
            
          case 'fake-jittering'
            for iJtr = 1 : statTestInfo.nJtr
                jtrd_spikeTrains = spikeTrains;
                fprintf('\n')
                [~,~, jtrdspt_gPLVs(iJtr),  ~, jtrdspt_couplingMatrix(:,:, iJtr) , jtrdspt_singularValues(:, iJtr)] = ...
                    gpla_core(jtrd_spikeTrains, lfpPhases, 0, iSV, [], plvNrmlzMethed, unwhitenOpr);
                fprintf('.')            
            end

            
          otherwise
            error('case is not considered')
        end
        
        % compute the p-value for the overal coupling
        pValue = nansum(jtrdspt_gPLVs >  gPLV) / statTestInfo.nJtr;
        gPLV_stats.nullDistribution = jtrdspt_gPLVs;
        nullHypoReject = pValue < statTestInfo.alphaValue;

        gPLV_stats.gPLV_testStat = gPLV;
        
        % compute the p-value for the individual spike-LFP couplings
        PLV_stats.pValue = ...
            nansum(abs(jtrdspt_couplingMatrix) ...
                   >  abs(repmat(couplingMatrix, 1, 1, statTestInfo.nJtr)), ...
                   3) ... % sum across the third dimention (repetetions)
            / statTestInfo.nJtr; % divide by number or repetion to get the p-value
        PLV_stats.nullDistribution = jtrdspt_couplingMatrix;
        PLV_stats.nullHypoReject = ...
            PLV_stats.pValue < statTestInfo.alphaValue * ones(size(PLV_stats.pValue));
        PLV_stats.nullDistribution_cPLV = jtrdspt_couplingMatrix;
        % SVspectrumStatsType = 'default';

        switch statTestInfo.SVspectrumStatsType

          case 'default'
            SVspectrum_stats.pValue = nansum(bsxfun(@gt, jtrdspt_singularValues, singularValues), 2) ...
                / statTestInfo.nJtr; % divide by number or repetion to get the p-value       
            SVspectrum_stats.nullDistribution = jtrdspt_singularValues;

          case 'cummulative'
            error('Need to be modfied for the new style of iSV input')
            % *** might need to be moved toward the end
            % nSV = numel(singularValues);
            % for iSV = 1 : nSV
            %     tmp_jtrdspt_singularValues = jtrdspt_singularValues(iSV:end, ...
            %                                                       :);
            %     SVspectrum_stats.pValue(iSV, 1) = ...
            %         nansum(sum(tmp_jtrdspt_singularValues .^ 2, 1) > sum(singularValues(iSV : end) .^ 2)) ...
            %         / statTestInfo.nJtr;

            % end
            % SVspectrum_stats.nullDistribution = jtrdspt_singularValues;
            
          case 'RMT-heuristic'
            % we replace the singular value spectrume with the largest ones
            % (repeating)
            jtrdspt_singularValues_onlyLargest = ...
                repmat(jtrdspt_singularValues(1,:), size(jtrdspt_singularValues, 1), 1);

            SVspectrum_stats.pValue = nansum(bsxfun(@gt, ...
                                                    jtrdspt_singularValues_onlyLargest, singularValues), 2) / statTestInfo.nJtr; 
            % divide by number or repetion to get the p-value       
           SVspectrum_stats.nullDistribution = jtrdspt_singularValues;

          otherwise
            error('no case for SV spec is considered')
        end
        
        varargout{4} = SVspectrum_stats;
      
      case 'RMT-based' % case on "testType"
        
        dimRatio = size(couplingMatrix, 1) / size(couplingMatrix, 2);
        lambda = (1 + dimRatio^.5) ^ 2;
        scaledSVspectrum = singularValues .^ 2 / size(couplingMatrix, 2);
        
        gPLV_stats.gPLV_testStat = scaledSVspectrum(1);
        nullHypoReject = (scaledSVspectrum(1) > lambda);
        gPLV_stats.nullHypoReject = nullHypoReject;
        
        SVspectrum_stats.nullHypoReject = (scaledSVspectrum > lambda);
        SVspectrum_stats.SV_testStat = scaledSVspectrum;
        
        pValue = NaN;
        PLV_stats = NaN;
        
        varargout{4} = SVspectrum_stats;
        
        % specify the significa
        % singular value spec
        
      otherwise % case on "testType"
        error('no case is considered')
    end

    % % (repeat a few stuff)
    % % check the data format
    % % data info
    % % specify the temporal window of analysis
    % % specify the units to be included in the analysis
    % % prepare the data

    %% Outputs 
    % pValue is already assigned        % 1     p-value
    varargout{1} = nullHypoReject;      % (2)   null hypothesis rejection
    varargout{2} = gPLV_stats;          % (3)   gPLV statistics (apart from the p-value)
    varargout{3} = PLV_stats;           % (4)   PLV statistics
                                        %                                   % (5)   SV statistics    
else
    pValue = NaN;
end

% if user asked for the normalized value the normalized value
% shoulod be returend
if flag_gPLVnrmlz == 1
    [~,~, gPLV] = ...
        gpla_core(spikeTrains, lfpPhases, flag_gPLVnrmlz, iSV, sameElecCheckInfo, ...
                  plvNrmlzMethed, unwhitenOpr);

    % gPLV = gPLV_out;
end