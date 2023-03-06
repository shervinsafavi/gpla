function [distances, varargout] = cmpt_arrayPWdistancesFromCenter(arrayMap)
% [distances, distancesTabel] = cmpt_arrayPWdistancesFromCenter(arrayMap)

ySize = size(arrayMap, 1);
xSize = size(arrayMap, 2);

refLocation = [round(ySize / 2) round(xSize / 2)];

% add warning if the size is such that one side will be slightly larger than the other one
% if (round(ySize / 2) || round(ySize / 2))

[distances, distancesTabel] = cmpt_arrayPWdistancesFromRef(arrayMap, refLocation)
varargout{1} = distancesTabel;    

% % refLocation : is the [xLocCenter, yLocCenter], where they are indices based on the array
% % ch =

% %     0.4057    0.1860    0.8261    0.9521    0.5594
% %     0.7999    0.8528    0.9386    0.8699    0.8209
% %     0.1017    0.3157    0.5366    0.0746    0.5141
% % center: [2, 3] ->  [yLocCenter, xLocCenter]

% % locDistanceUnit = vararargin{1};

% locDistanceUnit = 1;

% % distances = nan(max(arrayMap(:)), 1);

% [nyLoc, nxLoc] = size(arrayMap);
% nLoc = nxLoc * nyLoc;
% distances = nan(nLoc , 1);
% xLocs = nan(nLoc , 1);
% yLocs = nan(nLoc , 1);

% for ixLoc = 1 : nxLoc
%     for iyLoc = 1 : nyLoc
        
%         locID = arrayMap(iyLoc, ixLoc);
        
%         xLocs(locID) = ixLoc;
%         yLocs(locID) = iyLoc;        

%         distances(locID) = ...
%             norm([iyLoc ixLoc] - refLocation);
%         % compute the pair-wise distances
%         % distancesTabel(ixLoc, iyLoc) = ...
%         %     norm([integratedTabel(ixLoc, 1) integratedTabel(ixLoc, 2)]*locDistanceUnit ...
%         %     - [integratedTabel(iyLoc, 1) integratedTabel(iyLoc, 2)]*locDistanceUnit);
%         % PWdistances(iyLoc, ixLoc) = PWdistances(ixLoc, iyLoc);
%     end
% end

% distancesTabel = table(xLocs, yLocs, distances);

% distancesTabel.Properties.VariableNames = {'xLoc','yLoc','distanceFromRef'}
