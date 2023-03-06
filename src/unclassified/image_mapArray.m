function mappedVector = image_mapArray(imgVector, arrayMap, varargin)
% mappedVector = image_mapArray(imgVector, arrayMap, plotFlag)
% ------
% Input:
% 1      i1: 
% Output:
% 1      o1
%
% ------
% potential improvments:
% (1)dependes on the content of imgVector (e.g. complex, real val) use the
% proper image function. It has been used in old MUAimage
% (2) in new fig or not
% ------
% Code Info:
%   creation: 2014-08-17 by ShS -> shervin.safavi@gmail.com
%   modification:
%       $ 2016-10-23 support complex numbers for the image

if sum(~isnan(arrayMap(:))) ~= numel(imgVector)
    error('Provided array map is not matching the input')
end

%% initialization
nCh = numel(imgVector);
mappedVector = nan(size(arrayMap));

if nargin > 2
    plotFlag = varargin{1};   
end

% assigning the diffault values
if ~exist('plotFlag', 'var') plotFlag = 1; end

%% mapping the immage vector (imgVector)
for iCh = 1 : nCh
   [r, c] = find(arrayMap == iCh);
   mappedVector(r, c) = imgVector(iCh);
end

%% plot
if plotFlag
    if isreal(mappedVector)
        imagescnan(mappedVector, 'NanColor',[1 1 1]);
    else
        imaghuesc(angle(mappedVector),abs(mappedVector),[-pi pi], ...
                  [0 max(abs(mappedVector(:)))*.5], plotFlag);
%         circ_collbar([0 max(abs(mappedVector(:)))*.5])
    end
    axis square; box off; axis off
end