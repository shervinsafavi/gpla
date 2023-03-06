function cmap = circ_collbar_wThTick2(iScale, tickFontSize, tickLocs)
%CIRC_COLLBAR computes a circular color map and plots a circular colorbar
%cmap=CIRC_COLLBAR(iScale)
%
%  inputs
%
%	iScale - [min max] vector for the brightness range of the colorbar
%
%
%  outputs
%
%	cmap - returns the circular colormap (with hue from 0 to 1 and
%	saturation 0 and  brightness 1
%IMPORTANT: plots the colorbar only with no output arguments
%
%
% Author : Michel Besserve, MPI for Intelligent Systems, MPI for Biological Cybernetics, Tuebingen, GERMANY


if nargin<1
    iScale = [0.1 1];
    tickFontSize = 20;
    tickLocs = ...
        [5.5 0; ...
         -.5 5.5; ...
         -6.7 0; ...
         -.6 -5.5]
end
iVals=linspace(iScale(1),iScale(2),20);

theta=-pi:0.02:pi;huescale=linspace(0,1,length(theta));
x=cos(theta);
y=sin(theta);
d=1:length(theta);
Ring=linspace(1,5,length(iVals));
if nargout==0
    % figure()
    % figure('name', 'HSV Circular Colorbar', 'NumberTitle','off')
    for kring=1:(length(Ring)-1)
        
        cmap=(hsv2rgb([huescale' ones(length(huescale),1) iVals(kring)*ones(length(huescale),1)]));
        cMap3=permute(cmap,[3 1 2]);
%         patch([Ring(kring)*x Ring(kring+1)*fliplr(x)],[Ring(kring)*y Ring(kring+1)*fliplr(y)], [cMap3(1,:,1) cMap3(1,:,1)]);
    patch([Ring(kring)*x Ring(kring+1)*fliplr(x)],[Ring(kring)*y Ring(kring+1)*fliplr(y)], ...
        rand(1,630, 3));
        p = patch([Ring(kring)*x Ring(kring+1)*fliplr(x)],[Ring(kring)*y Ring(kring+1)*fliplr(y)], ...
            cat(2,cMap3,flipdim(cMap3,2)), 'EdgeColor', 'flat');
    end
    axis off
    axis equal
    set(gcf,'Renderer','zbuffer')
end
cmap=(hsv2rgb([huescale' ones(length(huescale),1) ones(length(huescale),1)]));

%%
% clf
% circ_collbar
% axis on

axis tight
xline(0);
yline(0);

%
text(tickLocs(1,1), tickLocs(1,2), '0', 'fontsize', tickFontSize, 'HorizontalAlignment', 'center')
text(tickLocs(2,1), tickLocs(2,2), '90', 'fontsize', tickFontSize, 'HorizontalAlignment', 'center')
text(tickLocs(3,1), tickLocs(3,2), '180', 'fontsize', tickFontSize, ...
     'HorizontalAlignment', 'center')
text(tickLocs(4,1), tickLocs(4,2), '270', 'fontsize', tickFontSize, ...
     'HorizontalAlignment', 'center')
