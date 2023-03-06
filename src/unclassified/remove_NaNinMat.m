function [mat_nanRemoved varargout] = remove_NaNinMat(mat, varargin)
% [mat_nanRemoved nans vals] = remove_NaNinMat(mat, NaNdistributionType, removeType)
%
%
% EXAMPLE:
%
%
% ------
% Input:
% 1     i1
%
% Output:
% 1     o1
%
% ------
% see also isnan
% ------
% potential improvments:
% (1) extending beyond 2D matrix
% (2) maybe write a new function which to the same thing with inf, empty
% guys collectivly
% ------
% Code Info:
%   creation: 2015-06-25 by ShS (shervin.safavi@gmail.com)
%   modification:
%       $ 201?-??-?? ?

%% handle optional inputs (varargin): NaNdistributionType, removeType
optionalVariables.NaNdistributionType = []; optionalVariables.removeType = [];
defaultValues{1} = 'sparse'; defaultValues{2} = 'both';
optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%%
% logicalInex_mat = isnan(mat);
% logicalInex_rows = sum(logicalInex_mat, 2);

[nR nC] = size(mat);

switch optionalVariables.NaNdistributionType
    case 'sparse'
        
        [tmp_nans.r tmp_nans.c] = find(isnan(mat));
        [tmp_vals.r tmp_vals.c] = find(~isnan(mat));
        
        nans.r = unique(tmp_nans.r);
        nans.c = unique(tmp_nans.c);
        
        vals.r = setdiff((1:nR), nans.r);
        vals.c = setdiff((1:nC), nans.c);
        
        % vals.r = unique(tmp_vals.r);
        % vals.c = unique(tmp_vals.c);
        
        switch optionalVariables.removeType
            
            case 'row'
                mat_nanRemoved = mat(vals.r, :);
            case 'col'
                mat_nanRemoved = mat(:, vals.c);
            case 'both'
                mat_nanRemoved = mat(vals.r, vals.c);
            otherwise
        end
        
    case 'fullRowOrCol'
        
        logicalIndex_mat = isnan(mat);
        
        % check which rows are fully NaN
        nNaNinRows = sum(logicalIndex_mat, 2);
        nans.r = find(nNaNinRows == nC);
        vals.r = setdiff((1:nR), nans.r);
        
        % check which columns are fully NaN
        nNaNinCols = sum(logicalIndex_mat, 1);
        nans.c = find(nNaNinCols == nR);
        vals.c = setdiff((1:nC), nans.c);
        
        switch optionalVariables.removeType
            case 'row'
                mat_nanRemoved = mat(vals.r, :);
            case 'col'
                mat_nanRemoved = mat(:, vals.c);
            case 'both'
                nans.bothRowAndCol = union(nans.r, nans.c);

                vals.r = setdiff((1:nR), nans.bothRowAndCol); 
                vals.c = setdiff((1:nC), nans.bothRowAndCol);    
                mat_nanRemoved = mat(vals.r, vals.c);
                
                
                % mat_nanRemoved = mat(:, vals.c);                                
                % % warning('If the all the row/col is NaN by removing in both dimentions you will have the empty matrix')                
            otherwise
                
        end
    otherwise
        
end


varargout{1} = nans;
varargout{2} = vals;