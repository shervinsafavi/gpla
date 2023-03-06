function relocateFigureStuff(axisHandle, relocationInfo)
% relocateFigureStuff(axisHandle, relocationInfo)
% relocationInfo = [new x-coordinate new y-coordinate width height]


pos = get(axisHandle, 'Position')

for k = 1 : numel(relocationInfo)
    if isnan(relocationInfo(k))
        updatedPosition(k) = pos(k);
    else
        updatedPosition(k) = relocationInfo(k);
    end
end

% posx = pos(1);
% posy = pos(2) + 0.01;
set(axisHandle, 'Position', updatedPosition)
