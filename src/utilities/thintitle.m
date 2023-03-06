function varargout = thintitle(varargin)

th = title(varargin{:});
set(th, 'FontWeight', 'normal');
varargout{1} = th;
