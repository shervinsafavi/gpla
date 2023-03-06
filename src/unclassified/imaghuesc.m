function ImgRes = imaghuesc(hMat, iMat, hueLim, intLim, varargin)
% ImgRes = imaghuesc(hMat, iMat, hueLim, intLim, plotFlag, colorbarFlag)
% IMAGHUESC hue valued matrix plot 
% plots the two input matrix as an image, the first matrix coding for hue
% and the second for intensity
%
% Inputs
%       hMat - hue value matrix
%       iMat - intensity (brightness) value matrix
%       hueLim - hue axis limits [min max]
%       intLim - intensity axis limits [min max]
%
% Outputs
%       ImgRes - RGB image output (lines x columns x 3)
%
%
% Author : Michel Besserve, 
% MPI for Intelligent Systems, MPI for Biological Cybernetics, Tuebingen, Germany

hMat(hMat>hueLim(2)) = hueLim(2);
hMat(hMat<hueLim(1)) = hueLim(1);
iMat(iMat>intLim(2)) = intLim(2);
iMat(iMat<intLim(1)) = intLim(1);

hMat=(hMat-hueLim(1)) / diff(hueLim);
iMat=(iMat-intLim(1)) / diff(intLim);

ImgRes = hsv2rgb(hMat, iMat*0+1, iMat);


%% plotting
if (nargin > 4), plotFlag = varargin{1}; end
if (exist('plotFlag', 'var') && plotFlag)
    % figure('name', ' - HSV Plot');
    imagescnan(ImgRes);
end

if (nargin > 5), colorbarFlag = varargin{2}; end
if (exist('colorbarFlag', 'var') && colorbarFlag)
   circ_collbar;
end

