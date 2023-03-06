function [singularLfpVecs, singularSpkVecs, singularValues] = fctrz_couplingMatrix(M)
% [singularLfpVecs, singularSpkVecs, singularValues] = fctrz_couplingMatrix(M)

%% SVD decompostion
if sum(isnan(M(:))) > 0 % means you had some units without any spike
    warning('You had units with no spikes (at least in some trials). These units were excluded from the analysis')
	[M_cleaned, ~, valIndex] = remove_NaNinMat(M, 'fullRowOrCol', 'col');
%     [~, c] = find(~isnan(M));
%     M_cleaned = M(:, c);
    [singularLfpVecs, D, V_cleaned] = svd(M_cleaned); 
    % singularValues = diag(D);
    singularSpkVecs = nan(size(M, 2));
    singularSpkVecs(valIndex.c,1:numel(valIndex.c)) = V_cleaned;
%     U(r, :)
else % normal procedure
    [singularLfpVecs, D, singularSpkVecs] = svd(M); 
    % singularValues = diag(D);
end

% if M is one-dimentional the diag will turn in into a matrix
if any(size(M) == 1)
    singularValues = D;
else
    singularValues = diag(D);
end