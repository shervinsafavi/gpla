%% MP distribution - with coupling
subplot2(nR,nC, 2, [3 4]);
ignit_gpla
load(fullfile(dirPath_gpla.src, 'visualizations', 'TPP_neuroIPSworkshop', 'multiDimSimulation04', 'dat', 'allVars.mat'))
% load /media/mpivision2/projects/shervin/research/gpla/src/visualizations/TPP_neuroIPSworkshop/multiDimSimulation04/dat/allVars.mat

for ic = 1 : ncase
    clear msd td_all

    caseName = caseNames{ic};
    tods = [];
    for iRel = 1 : nRel
        tmp = (svdOut.(caseName){icc, iRel}.singularValues) .^ 2 / signalParams.nUnit; % temp original data specturme
        tods = [tods; tmp];
    end

    histAx = linspace(0, max(tods)+1, 50);

    td = hist(tods, histAx) / numel(tods);
    % plot(histAx, td, 'r', 'linewidth', lw);
    bar(histAx, td, 'r')
    
    dimRatio = size(svdOut.(caseName){icc, iRel}.couplingMatrix, 1) / size(svdOut.(caseName){icc, iRel}.couplingMatrix, 2);
    lambda = (1 + dimRatio^.5) ^ 2;

    % vline(lambda, 'b', 'RMT-based threshold')

    % Ratio of matrix dimensions
    c = dimRatio;
    s = 1; % variance

    % Boundaries 
    a=(s^2)*(1-sqrt(c))^2;
    b=(s^2)*(1+sqrt(c))^2;

    dx = .05;
    svGrid = histAx;

    % Theoretical pdf
    ft=@(svGrid, a,b,c) (1./(2*pi* svGrid * c * s^(2))).*sqrt((b - svGrid).*(svGrid - a));
    F = real(ft(svGrid,a,b,c)); % for the values outside domain defined real value will be zero
    pdf = F/(sum(F));

    hold all
    plot(svGrid, pdf,'k','LineWidth', lw);
    % plot(svGrid, pdf,'ko','LineWidth',2);
    vline(lambda, 'b', 'RMT-based threshold')
    grid on

    % ylabel('Probability')
    xlabel('Scaled eigen value')
    % xlabel('Scaled singular value^2')
    % xlabel('Singular value^2 / num. unites')
end
set(gca, 'fontsize', fs)
legend('Empirical distribution', 'MP distribution')
