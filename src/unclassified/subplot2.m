% linearIndex = subplot2(ny, nx, py, px, anyOtherOption)
% 2D indexing of subplot nx, ny as before
% but px and py are subplot index coordinates 
% e.g. subplot(4,4,3) => subplot2(4,4,2,1);
% it is a simple enough function; I should mention that I saw it in Anton Sirota's lab
function [linearIndex, varargout] = subplot2(ny, nx, py, px, varargin)
handle = subplot(ny, nx, (py - 1) * nx + px, varargin{:});
% subplot(ny, nx, (py - 1) * nx + px, varargin);
linearIndex = (py - 1) * nx + px;

varargout{1} = handle;