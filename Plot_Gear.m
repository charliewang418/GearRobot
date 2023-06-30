function Plot_Gear(xyt, Epara)
%%
N = Epara.N;
Ng = Epara.Ng;
Ns = Epara.Ns;
Rg = Epara.Rg;
Rb = Epara.Rb;
Db = Epara.Db;
idx_start = Epara.idx_start;
idx_end = Epara.idx_end;

dtheta = 2 * pi / Ns;

x = zeros(N, 1);
y = zeros(N, 1);
for ng = 1:Ng
    xc = xyt(ng);
    yc = xyt(ng + Ng);
    theta_s = xyt(ng + 2 * Ng) + (0:Ns-1) * dtheta;
    idx1 = idx_start(ng);
    idx2 = idx_end(ng);
    x(idx1:idx2) = xc + Rg * cos(theta_s);
    y(idx1:idx2) = yc + Rg * sin(theta_s);
end
%%
figure; hold on; box on;

for ng = 1:Ng
    for i = idx_start(ng):idx_end(ng)
        rectangle('Position', [x(i) - Rb, y(i) - Rb, Db, Db],...
                  'Curvature', [1 1], 'EdgeColor', 'b', 'FaceColor', 'b');
    end
end

% plot([-1 2 * Ng], [0 0], 'k-', 'linewidth', 2)
for ns = round(-1 / Db):round(2 * Ng / Db)
    rectangle('Position', [(ns - 1) * Db - Rb, -Rb, Db, Db],...
              'Curvature', [1 1], 'EdgeColor', 'k', 'FaceColor', 'k');
end

axis equal