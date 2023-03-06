function [whitenDataMat, varargout] = whitenRed4(dataMat_raw,proportion)
% select pcas accouting for proportion*100 pct of the variance and whiten
% the signal by projecting on those (this also reduced dimension)
% inputs:
% * dataMat_raw: data matrix (channels/features x samples)
% * proportion: proportion of variance to account for, defaut .99
if nargin<2
    proportion = .99;
end
[tmpSig, whitenOpr, whitenInvOpr] = whitenRed2(reshape(dataMat_raw, size(dataMat_raw, 1), []),proportion);
nred = size(whitenOpr,1);
clear tmpSig
%whitenDataMat = reshape(tmpSig, size(tmpSig, 1), [], nTrial);    
nTrial = size(dataMat_raw,3);
whitenDataMat = zeros(nred,size(dataMat_raw,2),nTrial);
for ktrial = 1:nTrial
    [wdata, wOp, iwOp] = whitenRed2(dataMat_raw(:,:,ktrial),nred);
    whitenDataMat(:,:,ktrial) = wdata;
end
%[tmpSig, whitenOpr, whitenInvOpr] = whitenRed2(reshape(dataMat_raw, size(dataMat_raw, 1), []));
%whitenDataMat = reshape(tmpSig, size(tmpSig, 1), [], nTrial);    
whitenOpr = reshape(dataMat_raw,size(dataMat_raw,1),[])'\reshape(whitenDataMat,nred,[])';
whitenInvOpr = reshape(whitenDataMat,nred,[])'\reshape(dataMat_raw,size(dataMat_raw,1),[])';
varargout{1} = whitenOpr';
varargout{2} = whitenInvOpr';
