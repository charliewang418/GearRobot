function Fval = Force_Bumpy(xyt, Epara)
%%
Ng = Epara.Ng;
Ns = Epara.Ns;
N = Epara.N;
K = Epara.K;
F_grav = Epara.F_grav;
Rg = Epara.Rg;
% Rb = Epara.Rb;
Db = Epara.Db;
idx_start = Epara.idx_start;
idx_end = Epara.idx_end;

% Eval = 0;
dtheta = 2 * pi / Ns;

x = zeros(N, 1); % x coordinates of all bumps
y = zeros(N, 1); % y coordinates of all bumps
theta = zeros(N, 1); % orientation of each bump with respect to gear center
for ng = 1:Ng
    xc = xyt(ng);
    yc = xyt(ng + Ng);
    theta_s = xyt(ng + 2 * Ng) + (0:Ns-1) * dtheta;
    idx1 = idx_start(ng);
    idx2 = idx_end(ng);
    x(idx1:idx2) = xc + Rg * cos(theta_s);
    y(idx1:idx2) = yc + Rg * sin(theta_s);
    theta(idx1:idx2) = theta_s;
end

Fx = zeros(Ng, 1); % x force
Fy = zeros(Ng, 1); % y force
Ft = zeros(Ng, 1); % torque

% forces between two bumps that belong to two different gear particles
for ng = 1:Ng-1
    for ns = idx_start(ng):idx_end(ng)
        x1 = x(ns);
        y1 = y(ns);
        for mg = ng+1:Ng
            for ms = idx_start(mg):idx_end(mg)
                dx = x(ms) - x1;
                if abs(dx) < Db
                    dy = y(ms) - y1;
                    if abs(dy) < Db
                        d = sqrt(dx^2 + dy^2);
                        if d < Db
                            F = K * (Db / d - 1);
                            dFx = F * dx;
                            dFy = F * dy;
%                             Eval = Eval + 0.5 * K * (1 - d / Db)^2;  % cell-cell PE
                            Fx(ng) = Fx(ng) - dFx;
                            Fy(ng) = Fy(ng) - dFy;
                            Fx(mg) = Fx(mg) + dFx; % 3rd law
                            Fy(mg) = Fy(mg) + dFy;
                            Ft(ng) = Ft(ng) - Rg * (cos(theta(ns)) * dFy - sin(theta(ns)) * dFx);
                            Ft(mg) = Ft(mg) + Rg * (cos(theta(ms)) * dFy - sin(theta(ms)) * dFx);
                        end
                    end
                end
            end
        end
    end
end

% for ng = 1:Ng
%     for ns = idx_start(ng):idx_end(ng)
%         if y(ns) < Rb
%             dFy = K * (Rb - y(ns));
%             Fy(ng) = Fy(ng) + dFy;
%             Ft(ng) = Ft(ng) + Rg * cos(theta(ns)) * dFy;
%         end
%     end
% end

% forces between gear bumps and surface bumps
for ng = 1:Ng
    for ns = idx_start(ng):idx_end(ng)
        if y(ns) < Db
            dy = y(ns);
            dx1 = x(ns) - floor(x(ns) / Db) * Db;
            dx2 = dx1 - Db;
            d1 = sqrt(dx1^2 + dy^2);
            d2 = sqrt(dx2^2 + dy^2);
            if d1 < Db
                F1 = K * (Db / d1 - 1);
                dFx1 = F1 * dx1;
                dFy1 = F1 * dy;
                Fx(ng) = Fx(ng) + dFx1; % 3rd law
                Fy(ng) = Fy(ng) + dFy1;
                Ft(ng) = Ft(ng) + Rg * (cos(theta(ns)) * dFy1 - sin(theta(ns)) * dFx1);
            end
            if d2 < Db
                F2 = K * (Db / d2 - 1);
                dFx2 = F2 * dx2;
                dFy2 = F2 * dy;
                Fx(ng) = Fx(ng) + dFx2; % 3rd law
                Fy(ng) = Fy(ng) + dFy2;
                Ft(ng) = Ft(ng) + Rg * (cos(theta(ns)) * dFy2 - sin(theta(ns)) * dFx2);
            end
        end
    end
end

% gravity forces
Fy = Fy + F_grav * Ns;
for ng = 1:Ng
    Ft(ng) = Ft(ng) + F_grav * Rg * sum(cos(theta(idx_start(ng):idx_end(ng))));
end

Fval = [Fx; Fy; Ft];