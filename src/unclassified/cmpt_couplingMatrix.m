function couplingMatrix = cmpt_couplingMatrix(spkTrain, lfpPhase, varargin)
% add case in which get all the trials and sum them up and normalize it by
% the total number of spikes

%% Handle optional inputs (varargin):
optionalVariables.normalizationMethod     = [];     defaultValues{1} = 'nSpk';
optionalVariables.checkSameElecStuff_flag  ...
                                          = [];     defaultValues{2} = 0;
optionalVariables.statTestInfo            = [];     defaultValues{3} = [];

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

%%
switch optionalVariables.normalizationMethod
    
    case 'nSpk'
        %         F = exp(1i * lfpPhase);
        % compute the cross-covariance matrix
        tmp_couplingMatrix = F * spkTrain';
        
        % normalizing
        normalizingFactor = sum(spkTrain,2 ).^ -1;
        % take care of the case with zero spikes
        normalizingFactor(isinf(normalizingFactor))= nan;
        normalizingMatrix = repmat(normalizingFactor', nCh, 1);
        couplingMatrix = ( normalizingMatrix .* abs(tmp_couplingMatrix) ) .* exp(1i * angle(tmp_couplingMatrix));
        %         tmpMat = ( normalizingMatrix .* abs(tmp_couplingMatrix) ) .* exp(1i * angle(tmp_couplingMatrix));
        %         couplingMatrix = full(tmpMat);
        % suggestion: couplingMatrix =  normalizingMatrix .* tmp_couplingMatrix ;
        
        %% *** need to be modifed for a proper input
        if optionalVariables.checkSameElecStuff_flag == 1
            statTestInfo = optionalVariables.statTestInfo;
            jtrd_spkTrain = ...
                jitter_binarySparseSpikeTrain(spkTrain, statTestInfo.jitterWinWidth, statTestInfo.spkSF);
            jtrd_couplingMatrix = cmpt_couplingMatrix(full(jtrd_spkTrain), lfpPhase);
%             iJtr5 = 1
%             print('self-elec stuff is considered')
            clear jtrd_spkTrain
            
            if size(statTestInfo.spkU_lfpCh_cnvrtTabel, 2) == 2
                uchta = statTestInfo.spkU_lfpCh_cnvrtTabel;
                % tabel is given
                for iU = 1 : size(statTestInfo.spkU_lfpCh_cnvrtTabel,1)
                    couplingMatrix(uchta(iU,2), iU) = ...
                        jtrd_couplingMatrix(uchta(iU,2), iU);
                    %                             couplingMatrix(uchta(iU,2), uchta(iU,1)) = ...
                    %                                 jtrd_couplingMatrix(uchta(iU,2), uchta(iU,1));
                end
            elseif size(statTestInfo.spkU_lfpCh_cnvrtTabel, 2) == 1
                %     inices are just matching
            else
                %     error(()
            end
        end
        
%         if isstruct(optionalVariables.statTestInfo)
%             statTestInfo = optionalVariables.statTestInfo;
%             if (isfield(statTestInfo, 'stataliz_flag') && statTestInfo.stataliz_flag == 1)
%                 if isfield(statTestInfo, 'spkU_lfpCh_cnvrtTabel')
%                     jtrd_spkTrain = ...
%                         jitter_binarySparseSpikeTrain(spkTrain, statTestInfo.jitterWinWidth, statTestInfo.spkSF);
%                     jtrd_couplingMatrix = cmpt_couplingMatrix(full(jtrd_spkTrain), lfpPhase);
%                     clear jtrd_spkTrain
%                     
%                     if size(statTestInfo.spkU_lfpCh_cnvrtTabel, 2) == 2
%                         uchta = statTestInfo.spkU_lfpCh_cnvrtTabel;
%                         % tabel is given
%                         for iU = 1 : size(statTestInfo.spkU_lfpCh_cnvrtTabel,1)
%                             couplingMatrix(uchta(iU,2), iU) = ...
%                                 jtrd_couplingMatrix(uchta(iU,2), iU);
% %                             couplingMatrix(uchta(iU,2), uchta(iU,1)) = ...
% %                                 jtrd_couplingMatrix(uchta(iU,2), uchta(iU,1));
%                         end
%                     elseif size(statTestInfo.spkU_lfpCh_cnvrtTabel, 2) == 1
%                         %     inices are just matching
%                     else
%                         %     error(()
%                     end
%                 else
%                     %                error()necessary field is not provided
%                 end
%             end
%         elseif isempty(optionalVariables.statTestInfo)
%             % nothing need to be done
%         else
%             % error
%         end
        
        % %% *** need to be modifed for a proper input
        %         if isstruct(optionalVariables.sameElecStuff)
        %             sameElecStuff = optionalVariables.sameElecStuff;
        %             jtrd_spkTrain = ...
        %                 jitter_binarySparseSpikeTrain(spkTrain, jitterWinWidth, statTestInfo.spkSF);
        %             jtrd_couplingMatrix = cmpt_couplingMatrix(jtrd_spkTrain, lfpPhase);
        %
        %             if size(sameElecStuff.spkU_lfpCh_cnvrtTabel, 2) == 2
        %                 uchta = sameElecStuff.spkU_lfpCh_cnvrtTabel;
        %                 % tabel is given
        %                 for iU = 1 : size(sameElecStuff.spkU_lfpCh_cnvrtTabel,1)
        %                     couplingMatrix(uchta(iU,1), uchta(iU,2)) = ...
        %                         jtrd_couplingMatrix(uchta(iU,1), uchta(iU,2));
        %                 end
        %             elseif size(sameElecStuff.spkU_lfpCh_cnvrtTabel, 2) == 1
        %                 %     inices are just matching
        %             else
        %                 %     error(()
        %             end
        %         elseif optionalVariables.sameElecStuff == 0
        %             % nothing need to be done
        %         else
        %             % error
        %         end
        
        
        
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
        
end



% % % to normalize the phase locking values by the number of spikes in the
% % % spike train
% if normalizeSqrt
%     %normalizingFactor = (tmp_spikeTrains * ones(size(tmp_spikeTrains, 2), 1)).^ -.5;
%     normalizingFactor = sum(tmp_spikeTrains,2 ).^ -.5;
% else
%
% %     for kSpk = 1:nSpk
% %         for kTrial = 1:nTrial
% %             spktimes = find(tmp_spikeTrains(kSpk,:));
% %             randtimes = randomize_spikes(spktimes,nSamples);
% %             randtimes = randtimes(randtimes<=nSamples & randtimes>0);
% %             randM(:,kSpk,kTrial) = sum(tmp_F(:,ceil(randtimes)),2);
% %
% %         end
% %     end
% %     assumeZeroMean = 0;
% %     switch assumeZeroMean
% %         case 1
% %             normalizingMatrix = 1./sqrt(nanmean(abs(randM).^2,3));
% %         case 0
% %             normalizingMatrix = 1./nanstd(randM,[],3);
% %     end
% end
%
%
%
%
%
% % tmp_nrmlz_M = tmp_M;

% end

