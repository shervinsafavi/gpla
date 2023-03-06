function plot_eventRaster(eventTimes, varargin)
% plot_eventRaster(allEventTimes, plotSpecification, binSize, plotStyle)
% This function provide a [temporal] raster plot for event data, wherein
% each event will be indicated by a small vertical line. It can visulize
% single or multiple (= nDim) event time series. 
% A good example of such plot is spike raster plots for multiple trials
% ------
% Input:
% 1     eventTimes: 
%           D * T binary array 
%           D * T sparse matrix
%           D dimentional cell array
%       (no default value)
%       It should be your event data (e.g. spike train). It can be
%       binary array, or sparse matrix or a cell array.  
%           binary array: d*t array
%           sparse matrix: indicate the ones in binary array
%           cell array: cell{d, 1}, like conventional spike trains 
% (2)   plotSpecification: cell array
%       (default value is an empty cell)
%       For example you can indicate the color of discrete events small
%       vertical line by putting {'color','k'} as an optional input
% (3)   binSize: scalar
%       Bin size (in arbitrary unit) for the events (e.g. spikes) 
% 
% Output:
%
% ------
% potential improvments:
% (1) "plotSpecification" need to be passed wiser
% (2) automatic defining the x-axis limit
% (3) not always plot is a suitable tool for viz; sometimes line, sometimes
% imagesc some times line but don't know how this need to be figured out
% automatically
% ------
% Code Info:
%   creation: 2016-09-26 by ShS (shervin.safavi@gmail.com)
%   modification:
%       $ 2016-09-28 Use simple 'plot' syntax for the sake of speed
%       $ 2016-11-26 complete the section for cell array input
% ------
% see also rasterplot rasterplot rasterplot rasterplot_wColor image_disConData
% ------
% WORKING EXAMPLE:
% an example syntax (if it's needed)

%% handle optional variables
% if nargin > 1
%     nOptionalVar = 1;
%     plotSpecification = varargin(nOptionalVar+1 : end);
%     optionalInputValues = varargin(1 : nOptionalVar);
% else
%     optionalInputValues = varargin;
%     plotSpecification = {};
% end
optionalVariables.plotSpecification = [];       defaultValues{1} = {};
optionalVariables.binSize = [];                 defaultValues{2} = 1;
optionalVariables.plotStyle = [];               defaultValues{3} = 'line';
optionalVariables.opt_eventTime2binaryMat = []; defaultValues{4} = [];

 
optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%%
% thrForLightRaster = 20;

%% plotting
if iscell(eventTimes)    
    nDim = numel(eventTimes);
    for iUnit = 1 : nDim % loop on mumber of trials/cells/...
        tmp_eventTimes = eventTimes{iUnit};
        for iEvent = 1 : length(tmp_eventTimes)
            line([tmp_eventTimes(iEvent) tmp_eventTimes(iEvent)], [0 1]+iUnit-1, ...
            optionalVariables.plotSpecification{:});         
        end
    end
    set(gca,'YDir','reverse');
    % axis tight
    % ylim([0 nDim])
    
elseif issparse(eventTimes)
    % to be completed later
elseif isstruct(eventTimes) 
    % to be completed later
else % otherwise, [persumbaly] is a binary matrix
    eventTimes = logical(eventTimes);
    nDim = size(eventTimes, 1);
    [events_yCoordinate, events_xCoordinate] = find(eventTimes);
    
    % find has a wired behaviour with nDim = 1 and nDim > 1, try this two: 
    % [events_yCoordinate, events_xCoordinate] = find(rand(1, 1000) > 0.99); whos events_yCoordinate
    % [events_yCoordinate, events_xCoordinate] = find(rand(2, 1000) > 0.99); whos events_yCoordinate
    % you'll see the events_x/yCoordinate change from horizentally
    % organized to vertiacally organized!!
    if nDim == 1
        events_xCoordinate = events_xCoordinate';
        events_yCoordinate = events_yCoordinate';
    end
    nEvent = numel(events_yCoordinate);
    switch  optionalVariables.plotStyle
        
        case 'line'
        % if there is not many dim/channles/trials then use a line for each event 
        % i.e. if vertically plot is not very tight
        line([events_xCoordinate events_xCoordinate]' * optionalVariables.binSize, ...
            repmat([-.5; .5], 1, nEvent) + repmat(events_yCoordinate', 2, 1), ...
            optionalVariables.plotSpecification{:});
        set(gca,'YDir','reverse');
        % axis tight
        xlim([0 size(eventTimes, 2)*optionalVariables.binSize])
        
        case 'scatter'        
        % if there is many dim/channles/trials, then use a simpler graphical unit "."                 
        plot(events_xCoordinate * optionalVariables.binSize, events_yCoordinate, 'LineStyle','none', 'Marker','.', optionalVariables.plotSpecification{:});
        set(gca,'YDir','reverse');
        % axis tight
        xlim([0 size(eventTimes, 2)*optionalVariables.binSize])
%         ylim([0 nDim])
        
        case 'image'
            imagesc(eventTimes)
        otherwise
    end
end

end % end 

% function binaryMat = eventTime2binaryMat(eventTimes)
%     
% for iUnit = 1 : numel(eventTimes) % loop on mumber of trials/cells/...
%     binaryMat(iUnit, eventTimes{iUnit}) = 1;
% end
% end
%% drafts
% if iscell(eventTimes)
%     % to be completed later   
% elseif issparse(eventTimes)
%     % to be completed later    
% else % otherwise, [persumbaly] is a binary matrix
%     eventTimes = logical(eventTimes);
%     nDim = size(eventTimes, 1);
%     if nDim <= thrForLightRaster 
%         % if there is not many dim/channles/trials then use the 
%         % defult component for each event (line)
%         for iDim = 1 : nDim
%             eventTimes_singleDim = eventTimes(iDim, :);
%             eventTrain = find(eventTimes_singleDim == 1);
%             for iEvent = 1 : numel(eventTrain)
%                 line([eventTrain(iEvent) eventTrain(iEvent)], [0 1] + iDim-1, varargin{:});
%             end            
%         end
%     else % nDim > 20
%          % if there is many dim/channles/trials, then use a simpler graphical 
%          % defult component for each event (line)
%         hold all
%         for iDim = 1 : nDim
%             eventTimes_singleDim = eventTimes(iDim, :);
%             eventTrain = find(eventTimes_singleDim == 1);
%             nEvent = numel(eventTrain);
%             for iEvent = 1 : nEvent
%                 plot(eventTrain, iDim*ones(1, nEvent), '.', varargin{:});
%             end            
%         end         
%     end
%         
% end