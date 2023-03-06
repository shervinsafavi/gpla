function [couplingMatrix varargout] = tncmpt_couplingMatrix(spkTrain, lfpPhase, varargin)
% add case in which get all the trials and sum them up and normalize it by
% the total number of spikes

%% Handle optional inputs (varargin):
optionalVariables.normalizationMethod     = [];     defaultValues{1} = 'nSpk';
optionalVariables.sameElecCheckInfo  ...
    = [];     defaultValues{2} = [];

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%%
nCh = size(lfpPhase, 1);

%%
% *** we might need to run some check here:
% the phases are in range ([-pi pi])
% user want to use the envelope

if isreal(lfpPhase)
    F = exp(1i * lfpPhase);
else
    F = abs(lfpPhase) .* exp(1i * angle(lfpPhase));
end

varargout{1} = sum(spkTrain, 2);


%%
switch optionalVariables.normalizationMethod
    
  case {'nSpk', 'nSpk-square-root'}
    % compute the cross-covariance matrix
        tmp_couplingMatrix = F * spkTrain';

    % normalizing
        switch optionalVariables.normalizationMethod
          case {'nSpk'}
            normalizingFactor = sum(spkTrain,2 ).^ -1;
          case {'nSpk-square-root'}
            normalizingFactor = sum(spkTrain,2 ).^ -.5;
        end
        
    % take care of the case with zero spikes
        normalizingFactor(isinf(normalizingFactor))= nan;
        normalizingMatrix = repmat(normalizingFactor', nCh, 1);
        
        
        couplingMatrix = ( normalizingMatrix .* abs(tmp_couplingMatrix) ) .* exp(1i * angle(tmp_couplingMatrix));
    %         tmpMat = ( normalizingMatrix .* abs(tmp_couplingMatrix) ) .* exp(1i * angle(tmp_couplingMatrix));
    %         couplingMatrix = full(tmpMat);
    % suggestion: couplingMatrix =  normalizingMatrix .* tmp_couplingMatrix ;
        
    %% *** need to be modifed for a proper input
        
        if ~isempty(optionalVariables.sameElecCheckInfo)
        sameElecCheckInfo = optionalVariables.sameElecCheckInfo;
        jtrd_spkTrain = ...
            jitter_binarySparseSpikeTrain(spkTrain, sameElecCheckInfo.jitterWinWidth, sameElecCheckInfo.spkSF);
        jtrd_couplingMatrix = cmpt_couplingMatrix(full(jtrd_spkTrain), lfpPhase);
        %             iJtr5 = 1
        %             print('self-elec stuff is considered')
        clear jtrd_spkTrain
        
        if size(sameElecCheckInfo.spkU_lfpCh_cnvrtTabel, 2) == 2
            uchta = sameElecCheckInfo.spkU_lfpCh_cnvrtTabel;
            % tabel is given
            for iU = 1 : size(sameElecCheckInfo.spkU_lfpCh_cnvrtTabel,1)
                couplingMatrix(uchta(iU,2), iU) = ...
                    jtrd_couplingMatrix(uchta(iU,2), iU);
                %                             couplingMatrix(uchta(iU,2), uchta(iU,1)) = ...
                %                                 jtrd_couplingMatrix(uchta(iU,2), uchta(iU,1));
            end
        elseif size(sameElecCheckInfo.spkU_lfpCh_cnvrtTabel, 2) == 1
            %     inices are just matching
        else
            %     error(()
        end
        end

  case 'var1_theoretical'
    F = exp(1i * lfpPhase);
    %         % compute the cross-covariance matrix
    tmp_couplingMatrix = F * spkTrain';
    
    % normalizing
    normalizingFactor = sum(spkTrain,2 ).^ -.5;
    % take care of the case with zero spikes
    normalizingFactor(isinf(normalizingFactor))= nan;
    normalizingMatrix = repmat(normalizingFactor', nCh, 1);
    couplingMatrix = ( normalizingMatrix .* abs(tmp_couplingMatrix) ) .* exp(1i * angle(tmp_couplingMatrix));
    %         tmpMat = ( normalizingMatrix .* abs(tmp_couplingMatrix) ) .* exp(1i * angle(tmp_couplingMatrix));
    %         couplingMatrix = full(tmpMat);
    % suggestion: couplingMatrix =  normalizingMatrix .* tmp_couplingMatrix ;
    
  case 'var1_empirical'
    % *** note that in this case, all the trials need to be passed to
    % *** be careful after averaging the resukting gPLV the variance
    % will decrease again
  otherwise
    error('case is not considered')
end

