function jittered_binarySparseSpikeTrain ...
    = jitter_binarySparseSpikeTrain(binarySparseSpikeTrain, jitterWW, SF)
% jittered_binaryLogicalSpikeTrain = jitter_binaryLogicalSpikeTrain(binaryLogicalSpikeTrain, jitterWW, SF)

% *** check if the format is correct

%%
nSample    = size(binarySparseSpikeTrain, 2);      % number of samples

%% find the spike times
[rowInds, colInds] = find(binarySparseSpikeTrain == 1);
% because of strange ehaviour of find for 1D and 2D matrices, we need it!
if size(binarySparseSpikeTrain, 1) > 1
    spikeRowIndices = rowInds';
    spikeSampleInices = colInds';
else
    spikeRowIndices = rowInds;
    spikeSampleInices = colInds;
end

spikeTimes = spikeSampleInices / SF;

jittered_spikeTimes = ...
    (jitterWW * 2) * floor(spikeTimes / (jitterWW * 2)) + (jitterWW * 2) * rand(1 , length(spikeTimes));

%% output
jittered_spikeSampleInices = ceil(jittered_spikeTimes * SF);
jittered_binaryLogicalSpikeTrain = zeros(size(binarySparseSpikeTrain));

if ~isempty(jittered_spikeSampleInices)
    jittered_spikeSampleInices_refined = jittered_spikeSampleInices; 
    outRangeIdxs = find(jittered_spikeSampleInices < 0  | jittered_spikeSampleInices > nSample);    
    jittered_spikeSampleInices_refined(outRangeIdxs) = randi([1 nSample], 1, numel(outRangeIdxs),1);
%     jittered_SampleInices_refined = ...
%         jittered_spikeSampleInices(jittered_spikeSampleInices > 0 & jittered_spikeSampleInices <= nSample);
%     spikeRowIndices_refined = ...
%         spikeRowIndices(jittered_spikeSampleInices > 0 & jittered_spikeSampleInices <= nSample);

    
%     jittered_SampleInices_refined = ...
%         jittered_spikeSampleInices(jittered_spikeSampleInices > 0 & jittered_spikeSampleInices <= nSample);
%     spikeRowIndices_refined = ...
%         spikeRowIndices(jittered_spikeSampleInices > 0 & jittered_spikeSampleInices <= nSample);
%     numel(spikeSampleInices)
%     numel(spikeSampleInices) - numel(jittered_SampleInices_refined)

    % convert the 2D indices of spike sample indices to a linear index for
    % a matrix
%     tmpLinearIndForSpk = sub2ind(size(binarySparseSpikeTrain), spikeRowIndices_refined, jittered_SampleInices_refined);
    tmpLinearIndForSpk = sub2ind(size(binarySparseSpikeTrain), spikeRowIndices, jittered_spikeSampleInices_refined);
    jittered_binaryLogicalSpikeTrain(tmpLinearIndForSpk) = 1;    
end

jittered_binarySparseSpikeTrain = sparse(jittered_binaryLogicalSpikeTrain);

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
%         jittered_spikeIndices(jittered_spikeIndices > 0 & jittered_spikeIndices <= nSample);
% 
%     jittered_binaryLogicalSpikeTrain(jittered_spikeIndices_refined) = 1;        
% end
% 
% jittered_binarySparseSpikeTrain = sparse(jittered_binaryLogicalSpikeTrain);