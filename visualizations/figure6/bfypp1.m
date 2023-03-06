function bfypp1
% (rlabel)


thetaticks(0:45:315)

% [rmin, rmax] = get(gca, 'rlim')
rRange = get(gca, 'rlim')
tn = find_closestEvenNum1(rRange(2));

% tn = ceil(rRange(2));
% if mod(tn, 2) == 0, 
%     tn = tn; 
% else
%     tn = tn + 1;
% end

rticks([0 tn/2 tn]);
rlim([0 tn])
ax = gca;

rruler = ax.RAxis;
% rruler.Label.String = rlabel;

% ax.RAxisLocation = 67.5000;
ax.RAxisLocation = 22.5;

ax.GridAlpha = .2;
ax.LineWidth = 1.5;

% rruler = ax.RAxis;
% rruler.Label.String = rlabel;

% rruler.Label.String = 'Spikes / second';
