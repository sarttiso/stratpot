grad_vary = load('testdata/gradient_vary.csv');
trac_vary = load('testdata/trace_vary2.csv');

% generate toy grid
nint =  30;
x = linspace(-0.1,1.1,nint);
y = x;

[X,Y] = meshgrid(x,y);
X = reshape(X,numel(X),1);
Y = reshape(Y,numel(Y),1);

P = [X,Y];

% interpolate potential
Z_vary = stratpot2D(grad_vary, trac_vary, P);

% interpolate potential at bed contacts
[~, uidx] = unique(trac_vary(:,3));
Z_beds_vary = stratpot2D(grad_vary, trac_vary, trac_vary(uidx,1:2));

% define stratigraphy in potential coordinates
bed1_pot_vary = [0.8, Z_beds_vary(1);
                 0.8, Z_beds_vary(2)];
bed2_pot_vary = [0.8, Z_beds_vary(2);
                 0.8, Z_beds_vary(3)];
bed3_pot_vary = [0.8, Z_beds_vary(3);
                 0.8, Z_beds_vary(4)];

X = reshape(X,nint,nint);
Y = reshape(Y,nint,nint);
Zgrd_vary = reshape(Z_vary,nint,nint);

%% strata

% vary
% top bed
bedidx = (trac_vary(:,3) == 1) | (trac_vary(:,3) == 2);
bed1_vary = trac_vary(bedidx,1:2);
polidx = boundary(bed1_vary, 0.5);
poly1_vary = polyshape(bed1_vary(polidx,:));
% middle bed
bedidx = (trac_vary(:,3) == 2) | (trac_vary(:,3) == 3);
bed2_vary = trac_vary(bedidx,1:2);
polidx = boundary(bed2_vary, 0.5);
poly2_vary = polyshape(bed2_vary(polidx,:));
% bottom bed
bedidx = (trac_vary(:,3) == 3) | (trac_vary(:,3) == 4);
bed3_vary = trac_vary(bedidx,1:2);
polidx = boundary(bed3_vary, 0.7);
poly3_vary = polyshape(bed3_vary(polidx,:));

% stratigraphies
strat1 = [0.1, 0; 0.1, 1.0];
strat2 = [0.8, 0; 0.8, 1.0];

% intersections with units
bed1_strat1_vary = intersect(poly1_vary, strat1);
bed2_strat1_vary = intersect(poly2_vary, strat1);
bed3_strat1_vary = intersect(poly3_vary, strat1);

bed1_strat2_vary = intersect(poly1_vary, strat2);
bed2_strat2_vary = intersect(poly2_vary, strat2);
bed3_strat2_vary = intersect(poly3_vary, strat2);

%% plot 

col = linspecer(4);

% specify contour levels
lev = -0.1:0.1:1.1;
nlev = length(lev);

% specify contour colors
cmap = gray(20);
cmap = cmap(end-nlev+2:end, :);

contourf(X, Y, Zgrd_vary, lev)
hold on
colormap(cmap)

% first plot "beds"
% bed 1
plot(poly1_vary, 'FaceColor', col(1,:), 'FaceAlpha', 0.9)
% bed 2
plot(poly2_vary, 'FaceColor', col(2,:), 'FaceAlpha', 0.9)
% bed 3
plot(poly3_vary, 'FaceColor', col(3,:), 'FaceAlpha', 0.9)


% just lines
[C,h] = contour(X, Y, Zgrd_vary, lev, 'linecolor', 'k');
clabel(C, h, 'LabelSpacing', 100)

quiver(grad_vary(:,1), grad_vary(:,2), grad_vary(:,3), grad_vary(:,4), 0.3, 'k')
scatter(trac_vary(:,1),trac_vary(:,2), 50, col(trac_vary(:,3),:), ...
    'filled','MarkerEdgeColor', 'k')
% strats
% plot(strat1(:,1), strat1(:,2), '^-', 'color', 'k', ...
%     'markerfacecolor', 'k', 'linewidth', 2)
% plot(strat2(:,1), strat2(:,2), '^-', 'color', 'k', ...
%     'markerfacecolor', 'k', 'linewidth', 2)
xlabel('x', 'fontsize', 12)
ylabel('z', 'fontsize', 12)
title('Lateral Thickness Variability Example')

% label colorbar more appropriately
cbar = colorbar('ticks', lev);
cbar.Label.FontSize = 12;
cbar.Label.String = 'potential';
caxis([-0.1, 1.1])

axis equal

%% save fig

export_fig figures/stratpot_example_2D.png -m2 -transparent