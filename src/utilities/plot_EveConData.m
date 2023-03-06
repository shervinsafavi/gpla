function plot_EveConData(continuousData, eventData, varargin)
% plot_EveConData(continuousData, eventData, varargin)

%% handle optional variables
nOpVar = 0; % number of optional variable
nOpVar = nOpVar + 1; opVar.conPltSpecs   = []; defaultValues{nOpVar} = {};
nOpVar = nOpVar + 1; opVar.eventPltSpecs = []; defaultValues{nOpVar} = {};
nOpVar = nOpVar + 1; opVar.binSize       = []; defaultValues{nOpVar} = 1;
nOpVar = nOpVar + 1; opVar.rasterStyle   = []; defaultValues{nOpVar} = 'line';
 
opVar = handleVarargin(varargin, opVar, defaultValues);

%%
[nDim_cd, ~] = size(continuousData);
[nDim_ed, ~] = size(eventData);

if (nDim_cd == 1 && nDim_ed ~= 1)
    % make the continuousData possitive if needed
    if sum(continuousData < 0) > 0
        tmpShiftVal = -1 * min(continuousData); 
    else
        tmpShiftVal = 0;
    end
    continuousData = continuousData + tmpShiftVal;
    max_conData = max(continuousData);
    
    continuousData_scaled = -1 * continuousData * (nDim_ed / max_conData);
    continuousData_plot = continuousData_scaled - mean(continuousData_scaled) ...
        + nDim_ed / 2 + .5;
    plot(continuousData_plot , 'color','k', opVar.conPltSpecs{:});
    hold on
	plot_eventRaster(eventData, opVar.eventPltSpecs, opVar.binSize, opVar.rasterStyle)
    axis off
    box off
    
elseif (nDim_cd == 1 && nDim_ed == 1)
    
   %%
    % make the continuousData possitive if needed
    if sum(continuousData < 0) > 0
        tmpShiftVal = -1 * min(continuousData); 
    else
        tmpShiftVal = 0;
    end
    continuousData = continuousData + tmpShiftVal;
    max_conData = max(continuousData);
    
    continuousData_scaled = continuousData / max_conData;
%     continuousData_plot = continuousData_scaled;   
   %%
    plot(continuousData_scaled, 'LineWidth',3, 'Color','k', opVar.conPltSpecs{:})
%     axis off; box off
    
    % hold all
    tmpEventTime = find(eventData);
    for iEvent = 1 : length(tmpEventTime)
        line([tmpEventTime(iEvent) tmpEventTime(iEvent)], 1.2+[0 .2], ...
            'color','r', 'LineWidth',3, opVar.eventPltSpecs{:});
    end
%     vline(0,'b:', ['\beta = ' num2str(beta)])
%     ylim([0 .7])
    axis off; box off
else
    % to be completed
end