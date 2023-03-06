function jittered_binarySparseSpikeTrain = jitter_binarySparseSpikeTrainGroup(binarySparseSpikeTrain, jitterWW, SF)
%
nSample    = size(binarySparseSpikeTrain, 2);      % number of samples

nUnit = size(binarySparseSpikeTrain,1);
jittered_binaryLogicalSpikeTrain = zeros(nUnit,nSample);

winL = floor(jitterWW*2*SF);
nWin = floor(nSample/winL);
for kwin = 1:nWin
    randShift = randi(winL);
    jittered_binaryLogicalSpikeTrain(:,(1:winL)+(kwin-1)*winL) = circshift((binarySparseSpikeTrain(:,(1:winL)+(kwin-1)*winL)),randShift-1,2) ;
end
if size(binarySparseSpikeTrain,2)>winL*nWin
    randShift = randi(winL);
    jittered_binaryLogicalSpikeTrain(:,nWin*winL:end) = circshift(binarySparseSpikeTrain(:,nWin*winL:end),randShift-1,2) ;
end
jittered_binarySparseSpikeTrain = sparse(jittered_binaryLogicalSpikeTrain);
%jittered_binarySparseSpikeTrain = binarySparseSpikeTrain;