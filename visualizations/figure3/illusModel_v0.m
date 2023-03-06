% clf
d = 1;
ms = 15;
nlw = 3;
elw = 4;

ec = hsv(6);

plot(0, 1, 'ko', 'markersize',1.2*ms, 'linewidth',nlw+1)

hold all
for k = 1 : 3
    plot(k*d-d/2, 0, 'k^', 'markersize',ms, 'linewidth',nlw)
    line([0 k*d-d/2], [1 0], 'linewidth', rand*elw, 'color', ec(k, :));
    
    plot(-k*d+d/2, 0, 'k^', 'markersize',ms, 'linewidth',nlw)
    line([0 -k*d+d/2], [1 0], 'linewidth', rand*elw, 'color', ec(3+k, :));
end

% colorbar

axis off; box off