lvm = image_mapArray(angle(lfpVec.uate2010{ifreq}), utahMaps(iDataset).LFP.map, 0);
% svm = image_mapArray(angle(spkVec(:,iFreq)),utahMaps(iDataset).MU.map,0);

elecLoc = (1:10)*(4/10);

% subplot 122
% plot(elecLoc,nanmean(svm,1),'r')
% hold on
% for k = 1 : 10
%     plot(elecLoc,svm(k,:),'r.')
% end
% legend('Mean', 'Individual elec row')
% ylim([-pi pi])
% xlim([0 4.5])
% xlabel('Distance [mm]')
% % ylabel('Phase profile [rad]')

% subplot 121
% plot(elecLoc,nanmean(lvm,1),'b')
% hold on
% for k = 1 : 10
%     plot(elecLoc,lvm(k,:),'b.')
% end
% legend('Mean', 'Individual elec row')
% ylim([-pi pi])
% xlim([0 4.5])
% xlabel('Distance [mm]')
% ylabel('Phase profile [rad]')


%% fit
X = repmat(elecLoc, 10,1);
% mdl_lvm = fitlm(X(:), lvm(:))
mdl_lvm = fitlm(X(:), rad2deg(lvm(:)))

%% plot
% pos = get(sbh, 'Position')
% posx = pos(1);
% posy = pos(2) + 0.01;
% set(sbh,...
    % 'Position', [posx posy pos(3) pos(4)])

ph = plot(mdl_lvm);
box off

ylim(rad2deg([-pi pi]))
set(gca, 'ytick', [-180 0 180]);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.0355 0.035]);

set(gca, 'fontsize', vc.f5.fs_text)
set(ph(1), 'MarkerSize', 5); % marker size of LFP phase
set(ph(1), 'Marker', '+')
set(ph(1), 'color', vc.lfpv)

set(ph(2), 'linewidth', vc.f5.ppLW); % marker size of LFP phase
set(ph(2), 'color', [1 0 1])

set(ph(3), 'linewidth', vc.f5.ppLW); % marker size of LFP phase
set(ph(3), 'color', [1 0 1])

th = title('dummy');
delete(th)

legend off
lh = legend('LFP phases', 'Linear fit', 'Confidence bounds', ...
           'location','northwest');

displaceFigureStuff(lh, [-.002 .064 0 0]);
set(lh,...
    'box', 'on', ...
    'fontsize', vc.f5.sigTriangleLFS);


% pos = get(h,'Position')

% % top right
% posx = 0.6556;
% posy = 0.4640;
% legWidth = 0;

% set(h,...
%     'Position', [posx posy legWidth pos(4)], ...
%     'box', 'on', ...
%     'fontsize', vc.f5.sigTriangleLFS);


% title('LFP vec')
ylim(rad2deg([-pi pi]))
% xlim([0 4.5])
xlim([0 4])
xlabel('Horizental locations [mm]')
ylabel('Phase profile [deg]')

% *** to be done
% ytick
% tick dir
% location of the legend 
% location of the plot

