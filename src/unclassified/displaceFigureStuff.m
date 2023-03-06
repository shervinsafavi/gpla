function displaceFigureStuff(axisHandle, displacmentInfo)
% displaceFigureStuff(axisHandle, displacmentInfo)
% placmentInfo = [x-coordinate y-coordinate width height]


pos = get(axisHandle, 'Position')

for k = 1 : numel(displacmentInfo)
    if isnan(displacmentInfo(k))
        updatedPosition(k) = pos(k);
    else
        updatedPosition(k) = pos(k) + displacmentInfo(k);
    end
end

% posx = pos(1);
% posy = pos(2) + 0.01;
set(axisHandle, 'Position', updatedPosition)
