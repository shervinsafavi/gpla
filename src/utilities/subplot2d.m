function [linearIndex, varargout] = subplot2d(ny, nx, py, px, varargin)
% linearIndex = subplot2d(ny, nx, py, px, varargin)

if numel(py) == 1
    subplotHandle = subplot(ny, nx, (py - 1) * nx + px, varargin{:});
    % subplot(ny, nx, (py - 1) * nx + px, varargin);
    linearIndex = (py - 1) * nx + px;
elseif numel(py) > 1
    nRow = numel(py);
    for k = 1 : nRow
        tmp_py = py(k);
        tmpVec(k, :) = (tmp_py - 1) * nx + px;
        linearIndex = tmpVec(:);
        subplotHandle = subplot(ny, nx, linearIndex, varargin{:});
    end
end

varargout{1} = subplotHandle;
