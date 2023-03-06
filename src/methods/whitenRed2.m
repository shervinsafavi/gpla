function [whitenDataMat, varargout] = whitenRed2(dataMat_raw,proportion)
% select pcas accouting for proportion*100 pct of the variance and whiten
% the signal by projecting on those (this also reduced dimension)
% inputs:
% * dataMat_raw: data matrix (channels/features x samples)
% * proportion: proportion <1 of variance to account for, defaut .99, or
% number of components to keep (integer>=1)
if nargin<2
    proportion = .99;
end
%if proportion>1, error('proportion should be 0< ... <1'),end
% [whitenDataMat, varargout] = whiten(dataMat_raw)
N = size(dataMat_raw, 2);
dataMat = bsxfun(@minus, dataMat_raw, mean(dataMat_raw, 2));
covMat = dataMat * dataMat' / N;

[U, D] = eig(covMat);
% eigenvalues are not always sorted, here we sor
[d,ind] = sort(diag(D),'descend');
cumTrace = cumsum(d);
if proportion<1
    n = find(cumTrace>=proportion*cumTrace(end));
    n = n(1);
elseif isnan(proportion)
    n = size(diag(D), 1);
    % this option is added in order to be able to use the same
    % routines for the case where no reduction in the rank is asked
    % by user for the whitening.
    % in a sense this of witing the code, i.e. "proportion" can be both
    % ratio and the number component is somewhat ambigious. But
    % given that original way of coding was like this 
    % (see commit SHA 42a6ed15e2690d93d2bb8529d2d45283aa19792f)
else
    n = proportion; % was nasty bug
    % n = size(diag(D), 1);
end
ind = ind(1:n);
Ds = D(ind,ind);
Us = U(:,ind);

% whitenDataMat = (D .^ -.5) * U' * dataMat;
whitenDataMat = diag(diag(Ds) .^ -.5) * Us' * dataMat;

whitenOpr = diag(diag(Ds) .^ -.5) * Us';
whitenInvOpr = Us * diag(diag(Ds) .^ .5);
varargout{1} = whitenOpr;
varargout{2} = whitenInvOpr;
varargout{3} = mean(dataMat_raw, 2);
