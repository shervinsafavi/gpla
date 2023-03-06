function jittered_binarySparseSpikeTrain = jitter_binarySparseSpikeTrainISI(binarySparseSpikeTrain, jitterWW, SF)
%%
nSample    = size(binarySparseSpikeTrain, 2);      % number of samples

nUnit = size(binarySparseSpikeTrain,1);
jittered_binaryLogicalSpikeTrain = zeros(nUnit,nSample);
for kunit = 1: nUnit
    %% find the spike times
    [spikeSampleInices] = find(binarySparseSpikeTrain(kunit,:) == 1);
    spikeTimes = spikeSampleInices / SF;

    % shuffle by peruting ISIs
    jitterWinIdx = floor(spikeTimes / (jitterWW * 2));
    for kwin = unique(jitterWinIdx)
        ISI = diff(spikeTimes(jitterWinIdx==kwin));
        % generate new times, relative to an initial time to be determined,
        % by permuting ISIs
        relativeNewTimes = [0,cumsum(ISI(randperm(length(ISI))))];
        % range we have to set initial spike while staying in the window
        initSpikeRange = 2*jitterWW-relativeNewTimes(end);
        spikeTimes(jitterWinIdx==kwin) = relativeNewTimes+initSpikeRange*rand(1)+(jitterWW *2)*kwin;
    end
    % add a jitter at the scale of binning period to avoid potential
    % discretization effects
    spikeTimes = spikeTimes+1/SF*rand(size(spikeTimes,1),size(spikeTimes,2));
    %% output
    jittered_spikeSampleInices = ceil(spikeTimes * SF);
    

    if ~isempty(jittered_spikeSampleInices)
        jittered_spikeSampleInices_refined = jittered_spikeSampleInices; 
        outRangeIdxs = find(jittered_spikeSampleInices < 0  | jittered_spikeSampleInices > nSample);    
        jittered_spikeSampleInices_refined(outRangeIdxs) = randi([1 nSample], 1, numel(outRangeIdxs),1);

        %tmpLinearIndForSpk = sub2ind(size(binarySparseSpikeTrain), spikeRowIndices, jittered_spikeSampleInices_refined);
        jittered_binaryLogicalSpikeTrain(kunit, jittered_spikeSampleInices_refined) = 1;    
    end
end
jittered_binarySparseSpikeTrain = sparse(jittered_binaryLogicalSpikeTrain);
