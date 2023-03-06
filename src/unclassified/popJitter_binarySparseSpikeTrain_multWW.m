function jittered_binarySparseSpikeTrain ...
    = popJitter_binarySparseSpikeTrain_multWW(binarySparseSpikeTrain, jitterWWrange, SF)
%   popJitter_binarySparseSpikeTrain_multWW
% *** check if the format is correct
% *** change the name to syncJitter...
%%
nSample    = size(binarySparseSpikeTrain, 2);      % number of samples

jitterWW = jitterWWrange(1) + diff(jitterWWrange) .* rand;

%% find the spike times
[rowInds, colInds] = find(binarySparseSpikeTrain == 1);
% because of strange ehaviour of find for 1D and 2D matrices, we need it!
if size(binarySparseSpikeTrain, 1) > 1
    spikeRowIndices = rowInds';
    spikeSampleInices = colInds';
else
% *** if this is the case doesnt make sense to used this function at all
    spikeRowIndices = rowInds;
    spikeSampleInices = colInds;
end

[popSpikeSampleInices,~, indexW_popSpikeSI]=  unique(spikeSampleInices); 
% indexW_popSpikeSI = indexWith_popSpikeSampleIices

popSpikeTimes = unique(popSpikeSampleInices) / SF;

jittered_popSpikeTimes = ...
    (jitterWW * 2) * floor(popSpikeTimes / (jitterWW * 2)) + (jitterWW * 2) * rand(1 , length(popSpikeTimes));

%% output
jittered_popSpikeSampleInices = round(jittered_popSpikeTimes * SF);
jittered_spikeSampleInices = jittered_popSpikeSampleInices(indexW_popSpikeSI);

jittered_binaryLogicalSpikeTrain = zeros(size(binarySparseSpikeTrain));

if ~isempty(jittered_spikeSampleInices)
    jittered_SampleInices_refined = ...
        jittered_spikeSampleInices(jittered_spikeSampleInices > 0 & jittered_spikeSampleInices <= nSample);
    spikeRowIndices_refined = ...
        spikeRowIndices(jittered_spikeSampleInices > 0 & jittered_spikeSampleInices <= nSample);
    % convert the 2D indices of spike sample indices to a linear index for
    % a matrix
    tmpLinearIndForSpk = sub2ind(size(binarySparseSpikeTrain), spikeRowIndices_refined, jittered_SampleInices_refined);
    jittered_binaryLogicalSpikeTrain(tmpLinearIndForSpk) = 1;    
end

jittered_binarySparseSpikeTrain = sparse(jittered_binaryLogicalSpikeTrain);

% % popSpikeTimes 
% popSpikeSampleInices
% full(binarySparseSpikeTrain)
% % jittered_popSpikeTirmes
% jittered_popSpikeSampleInices
% full(jittered_binarySparseSpikeTrain)
% 
% %% find the spike times
% spikeIndices = find(binarySparseSpikeTrain == 1);
% spikeTimes = spikeIndices / SF;
% 
% jittered_spikeTimes = ...
%     (jitterWW * 2) * floor(spikeTimes / (jitterWW * 2)) + (jitterWW * 2) * rand(1 , length(spikeTimes));
% 
% %% output
% jittered_spikeIndices = round(jittered_spikeTimes * SF);
% jittered_binaryLogicalSpikeTrain = zeros(size(binarySparseSpikeTrain));
% if ~isempty(jittered_spikeIndices)
%     jittered_spikeIndices_refined = ...
%         jittered_spikeIndices(jittered_spikeIndices > 0 & jittered_spikeIndices <= nSampl
%     jittered_binaryLogicalSpikeTrain(jittered_spikeIndices_refined) = 1;        
% end
% 
% jittered_binarySparseSpikeTrain = sparse(jittered_binaryLogicalSpikeTrain);
