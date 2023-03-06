% clf
% [~, sbh] = subplot2(nR,nC, nred+1, [3 4]);

% nSelFreq = 4;
nSelFreq = 3;

allColors = cool(nSelFreq);

clear ph

for kFreq = 1 : nSelFreq
    % scatter(abs(spkVec.uate2010{kFreq}), rad2deg(angle(spkVec.uate2010{kFreq})), 'filled', 'MarkerFaceColor', allColors(kFreq, :));
    xVals = abs(spkVec.uate2010{kFreq});
    scatter(xVals, rad2deg(angle(spkVec.uate2010{kFreq})), 3, 'MarkerFaceColor', allColors(kFreq, :));
    hold on
    
    mdl_lvm = fitlm(abs(spkVec.uate2010{kFreq}(:)), rad2deg(angle(spkVec.uate2010{kFreq}(:))));
   
    yFit = xVals * mdl_lvm.Coefficients.Estimate(2) + mdl_lvm.Coefficients.Estimate(1);
    ph(kFreq) = plot(xVals, yFit,'color', [1 0 0],'Linewidth',2, 'color', allColors(kFreq, :))
    
    % plot(x,x*res.Coefficients.Estimate(2)+res.Coefficients.Estimate(1),'color',h1.Color,'Linewidth',2)
    % title(freqBandsLab{kFreq})
    % if (any(kFreq == [1 5])), ylabel('Phase [deg]'); end 
    % if (any(kFreq == (5:nFreq))), xlabel('Modulus [a.u.]'); end 
end

ylim(100 * [-1 1])
set(gca, 'ytick', [-100 0 100]);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.0355 0.035]);

set(gca, 'fontsize', vc.f5.fs_text)

% lh = legend(ph, 'LFP phases', 'Linear fit', 'Confidence bounds', 'sdf', ...
%            'location','northwest');

lh = legend(ph, freqBandsLabHz(1:nSelFreq), ...
           'location','northwest');



% displaceFigureStuff(lh, [-.002 .064 0 0]);
displaceFigureStuff(lh, [-.002 .058 0 0]);
set(lh,...
    'box', 'off', ...
    'fontsize', vc.f5.sigTriangleLFS);

xlabel('Modulus [a.u.]')
ylabel('Phase [deg]')


%% storing model fitting 
% *** save the mdl_lvm 
clear mdl_lvm

for kFreq = 1 : nSelFreq

    xVals = abs(spkVec.uate2010{kFreq});
    mdl_lvm{kFreq} = fitlm(abs(spkVec.uate2010{kFreq}(:)), rad2deg(angle(spkVec.uate2010{kFreq}(:))));

end

%%
mdl_lvm{2}.Coefficients.pValue(2) < 10e-6
