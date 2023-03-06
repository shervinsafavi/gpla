% x = 1 : 100;
% y = 

%% v1
% clf

linewidth = 4.5; % desired linewidth
tau = .3;

x = 0:.1:10; % x data. Assumed to be increasing
y = 1 + exp(-tau*abs(x)); % y data
N = 100; % number of colors. Assumed to be greater than size of x
cmap = cool(N); % colormap, with N colors
xi = x(1)+linspace(0,1,N+1)*x(end); % interpolated x values
yi = interp1(x,y,xi); % interpolated y values
hold on

% N = 50;


for n = 1:N
    plot(xi([n n+1]), yi([n n+1]), 'color', cmap(n,:), 'linewidth', linewidth);
end
hold on 
for n = 1:N
    plot(-xi([n n+1]), yi([n n+1]), 'color', cmap(n,:), 'linewidth', linewidth);
end

% x = 0:.1:10; % x data. Assumed to be increasing
y = exp(-tau*abs(x)); % y data
N = 100; % number of colors. Assumed to be greater than size of x
cmap = cool(N); % colormap, with N colors
% linewidth = 4.5; % desired linewidth
xi = x(1)+linspace(0,1,N+1)*x(end); % interpolated x values
yi = interp1(x,y,xi); % interpolated y values
hold on

% N = 50;

cmapR = flip(cmap);

for n = N:-1:1
    plot(xi([n n+1]), yi([n n+1]), 'color', cmapR(n,:), 'linewidth', linewidth);
end
for n = N:-1:1
    plot(-xi([n n+1]), yi([n n+1]), 'color', cmapR(n,:), 'linewidth', linewidth);
end

% colorbar
cbh = colorbar
caxis([0 180])
th = title(cbh, 'Phase lag');
displaceFigureStuff(th, [NaN -285 NaN])
colormap(gca, 'cool');
set(cbh, 'tickdir','out')
set(cbh, 'xtick', [0 90 180]);
% set(cbh, 'xtick', []);

% ylim([-.1 .7])
axis off

vc.f13.schfs = 7;
vc.f13.schfsArrow = 8;

xline(0, 'k:', 'linewidth', 2.5)
text(0, -.1, 'Center', 'HorizontalAlignment', 'center', 'FontWeight','normal', 'fontsize', vc.f13.schfs)
text(-14.5, 1.8, 'Weak recurrence', 'HorizontalAlignment', 'center', 'FontWeight','normal', 'fontsize', vc.f13.schfs)
text(-14.5, .8, 'Strong recurrence', 'HorizontalAlignment', 'center', 'FontWeight','normal', 'fontsize', vc.f13.schfs)

set(cbh, 'fontsize', vc.f13.schfs)


% x = [.15 .25];
% y = .25 * ones(1,2);
% annotation('textarrow',x,y,'String',' Gradient  ','FontSize',vc.f13.schfsArrow,'Linewidth',2)
% x = [.75 .65];
% annotation('textarrow',x,y,'String','','FontSize',vc.f13.schfsArrow,'Linewidth',2)

% x = flip([.15 .25]);
% y = .67 * ones(1,2);
% annotation('textarrow',x,y,'String','','FontSize',vc.f13.schfsArrow,'Linewidth',2)
% x = flip([.75 .65]);
% annotation('textarrow',x,y,'String','','FontSize',vc.f13.schfsArrow,'Linewidth',2)





