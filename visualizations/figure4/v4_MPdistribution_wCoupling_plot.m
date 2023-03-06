%% MP distribution - with coupling
subplot2(nR,nC, 2, [3 4]);
icc = 1; 

%%
run v3_commonStuff.m
load(fullfile('dat', 'v4_MPdistribution_wCoupling.mat'))

%%
clear msd td_all

tods = [];
for iRel = 1 : nRel
    tmp = (svdOut.(caseName){icc, iRel}.singularValues) .^ 2 / signalParams.nUnit; % temp original data specturme
    tods = [tods tmp];
end

histAx = linspace(0, max(tods)+1, 50);

%%
td = hist(tods, histAx) / numel(tods);
% plot(histAx, td, 'r', 'linewidth', lw);
b = bar(histAx, td, 'facecolor',.5*ones(1,3), 'EdgeColor',.4*ones(1,3), ...
        'BarWidth',.9)

dimRatio = size(svdOut.(caseName){icc, iRel}.couplingMatrix, 1) / size(svdOut.(caseName){icc, iRel}.couplingMatrix, 2);
lambda = (1 + dimRatio^.5) ^ 2;


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
plot(svGrid, pdf,'r','LineWidth', lw);
vline(lambda, 'b', 'Theoretical threshold')
axis tight

% plot(svGrid, pdf,'ko','LineWidth',2);

%%
grid on

% ylabel('Probability')
xlabel('Scaled eigen value')
% xlabel('Scaled singular value^2')
% set(gca, 'fontsize', fs)
% legend('Empirical distribution', 'MP distribution', 'location','northeastoutside')
legend('Empirical distribution', 'MP distribution')
