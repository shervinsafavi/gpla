function plotScaleBar1(X, Y, scaleLabel_x, scaleLabel_y, ...
                       xScale_lineLength, yScale_lineLength, ...
                       ratio_forYoffSet, ratio_forXoffSet, ... 
                       ratio_forTextYloc, ratio_forTextXloc, ...
                       LineWidth,  fs, varargin)
% plotScaleBar1
% should be organized row-wise

% unit.X = 'mm';
% unit.Y = '';
% xRange = 1;
% yRange = 0.02;
% LineWidth = 5;

% plot(data(1, :), data(2, :), varargin{:})
axis off; box off
% % assuming your plot is already there, so what ever is there should be kept
% % and the scale line should be plotted on top of it, so:
hold all
% minmax can be used
minX = min(X);
maxX = max(X)
xRange = max(X) - minX;

minY = min(Y);
yRange = max(Y) - minY ;


% X and Y ranges
% iDim = 1;
% minX = min(data(iDim, :));
% xRange = max(data(iDim, :)) - minX;
if ~isnan(scaleLabel_x)
    % iDim = 2;
    % minY = min(data(iDim, :));
    % yRange = max(data(iDim, :)) - minY ;

    % X scale
    xScale_lineLength_fd = .2 * xRange; % fd stands for 'from data'
    xScale_yOffset = ratio_forYoffSet * yRange;

    % X_xScaleLocation = [minX minX+xScale_lineLength];
    % Y_xScaleLocation = (minY -  xScale_yOffset) * [1 1];
    
    X_xScaleLocation = [maxX-xScale_lineLength maxX];
    Y_xScaleLocation = (minY -  xScale_yOffset) * [1 1];

    line(X_xScaleLocation, Y_xScaleLocation, 'LineWidth',LineWidth, 'color', 'k')
    % textUndeScaleLine = [num2str(xScale_lineLength), unit.X];
    textUndeScaleLine = scaleLabel_x;
    textLoc = max(X) - xScale_lineLength/2;
    % minY-2*xScale_yOffset
    % 5*xScale_yOffset = 5 *
    xScaleText_yOffset = ratio_forTextYloc * yRange;
    text(textLoc, minY - xScale_yOffset - xScaleText_yOffset, textUndeScaleLine, ...
         'FontSize', fs, 'FontName','SansSerif', ...%'FontWeight','bold', ...
         'color', 'k', 'HorizontalAlignment','center', ...
         'interpreter', 'none');

end

if ~isnan(scaleLabel_y)
    % iDim = 2;
    % minY = min(data(iDim, :));
    % yRange = max(data(iDim, :)) - minY ;

    % Y scale
    yScale_lineLength_fd = .2 * yRange; % fd stands for 'from data'
    yScale_xOffset = ratio_forXoffSet * xRange;
    % yScale_xOffset = ratio_forXoffSet * yRange;
    
    % X_xScaleLocation = [minX minX+xScale_lineLength];
    % Y_xScaleLocation = (minY -  xScale_yOffset) * [1 1];
    
    Y_yScaleLocation = [minY minY+yScale_lineLength];
    X_yScaleLocation = (maxX +  yScale_xOffset) * [1 1];

    line(X_yScaleLocation, Y_yScaleLocation, 'LineWidth',LineWidth, 'color', 'k')
    % textUndeScaleLine = [num2str(xScale_lineLength), unit.X];
    textUndeScaleLine = scaleLabel_y;
    % textLoc = min(Y) + yScale_lineLength/2;
    textLoc = maxX +  yScale_xOffset / 2;
    % minY-2*xScale_yOffset
    % 5*xScale_yOffset = 5 *
    yScaleText_xOffset = ratio_forTextXloc * xRange;
    % yScaleText_xOffset = ratio_forTextXloc * yRange;
    % text(textLoc, minY - yScale_xOffset - yScaleText_xOffset, textUndeScaleLine, ...
    %      'FontSize', fs, 'FontName','SansSerif', ...%'FontWeight','bold', ...
    %      'color', 'k', 'VerticalAlignment','middle', ...
    %      'interpreter', 'none');
    text(textLoc + yScaleText_xOffset + yScale_xOffset, minY + yScale_lineLength / 2 , textUndeScaleLine, ...
         'FontSize', fs, 'FontName','SansSerif', ...%'FontWeight','bold', ...
         'color', 'k', 'VerticalAlignment','middle', ...
         'interpreter', 'none');
    

end

% if ~isnan(scaleLabel_y)
%     error('not ready yet')
    
%     % % minY = min(Y);
%     % % yRange = max(Y) - minY ;

%     % % y scale
%     % yScale_lineLength_fd = .2 * yRange; % fd stands for 'from data'
%     % yScale_xOffset = .05 * xRange;

%     % X_yScaleLocation = (minX -  yScale_xOffset) * [1 1];
%     % Y_yScaleLocation = [minY minY+yScale_lineLength];

%     % line(X_yScaleLocation, Y_yScaleLocation, 'LineWidth',LineWidth, 'color', 'k')

%     % % textUndeScaleLine = [num2str(yScale_lineLength), unit.Y];
%     % textUndeScaleLine = scaleLabel_y;
%     % warning('check the line below to fix it similar to X one')
%     % text(minX-3*yScale_xOffset, minY+yScale_lineLength/2, textUndeScaleLine, ...
%     %      'FontSize', fs, 'FontName','SansSerif', ... %'FontWeight','bold', ...
%     %      'color', 'k', 'HorizontalAlignment','center', ...
%     %      'interpreter', 'none');

% end

