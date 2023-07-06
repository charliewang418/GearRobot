function GearAssembly(Ng, Ns, pid)
%% Ng: number of gear particles
%% Ns: number of bumps per gear particle
%% pid: seed for random number generator
%% Initialization
N = Ng * Ns; % total number of bumps
K = 1; % spring constant for all contacts
a_cut_c = 1.4;
Ka_ratio = 0.2;
F_grav = -0.0001; % gravity force on each bump
dt_od = 0.01; % time size for damped molecular dynamics (MD) to find force balanced state
Nt_od = 1e6; % max number of time steps for damped MD
Fthresh = 1e-13; % force balance threshold

idx_start = (0:Ng-1)' * Ns + 1; % starting index for bumps in a gear particle
idx_end = (1:Ng)' * Ns; % ending index for bumps in a gear particle

Rg = 0.5; % radius of the gear circle
Rb = Rg * sin(pi / Ns); % radius of the bump
Db = 2 * Rb; % diameter of the bump

Epara = struct('Rg', Rg, 'Rb', Rb, 'Db', Db, 'Ng', Ng, 'Ns', Ns, 'N', N,...
               'idx_start', idx_start, 'idx_end', idx_end,...
               'F_grav', F_grav, 'K', K, 'Ka_ratio', Ka_ratio, 'a_cut_c', a_cut_c);

xyt_init = GearInitialize(pid, Epara); % initialize configuration
% energy minimize the configuration so that particles are in force and torque balance
xyt_init = Damped_MD_Minimization(xyt_init, Epara, Fthresh, dt_od, Nt_od); 
%% Motion Test -- apply a constant torque, with damping
dt = 0.01; % time size for constant energy MD
Nt = 4e5; % total number of time steps for MD
dt_half = dt / 2;
dt_sq_half = 0.5 * dt * dt;
B = 0.05; % damping coefficient, set to 0 for constant energy MD
B_denom = 1 + B * dt_half;

T_drive = 0.003;
if Ng == 1
    T_ext = -T_drive;
else
    if mod(Ng, 2) == 0
        T_ext = T_drive * [-ones(Ng / 2, 1); ones(Ng / 2, 1)]; % applied constant torque on all particles
    else
        T_ext = T_drive * [-ones((Ng - 1) / 2, 1); 0; ones((Ng - 1) / 2, 1)];
    end
end
ext_idx1 = 2 * Ng + 1; % starting index for torque in the force&torque array
ext_idx2 = 3 * Ng; % ending index for torque in the force&torque array

xyt = xyt_init;
t = 0;

nt_record = 100; % store particle positions every nt_record time step
xyt_record = zeros(3 * Ng, Nt / nt_record + 1); % array to store particle positions during MD simulation
xyt_record(:, 1) = xyt;

Vel = zeros(3 * Ng, 1); % translation and rotation velocity
Acc = Force_Bumpy(xyt, Epara); % force and torque calculation
Acc(ext_idx1:ext_idx2) = Acc(ext_idx1:ext_idx2) + T_ext; % add applied external torque
Acc_old = Acc;

% MD using Verlet method
for nt = 1:Nt
    % update position by a full step
    xyt = xyt + Vel * dt + Acc_old * dt_sq_half;

    if mod(nt, nt_record) == 0
        xyt_record(:, nt / nt_record + 1) = xyt;
    end
    if mod(nt, 100000) == 0
        fprintf('nt: %d\n', nt);
    end

    % update force by a full step
    t = t + dt;
    Acc = Force_Bumpy(xyt, Epara);
    Acc(ext_idx1:ext_idx2) = Acc(ext_idx1:ext_idx2) + T_ext;
    % correct force with damping term
    Acc = (Acc - B * (Vel + Acc_old * dt_half)) / B_denom;
    % update velocity by a full step
    Vel = Vel + (Acc + Acc_old) * dt_half;
    Acc_old = Acc;
end

%% Motion Test -- apply oscillating torques, with damping
dt = 0.01; % time size for constant energy MD
Nt = 1e6; % total number of time steps for MD
dt_half = dt / 2;
dt_sq_half = 0.5 * dt * dt;
B = 0.05; % damping coefficient, set to 0 for constant energy MD
B_denom = 1 + B * dt_half;

T_drive = 0.006;
f_drive = 0.1;
if Ng == 1
    theta_drive_init = 3 * pi / 2;
else
%     if mod(Ng, 2) == 0
%         T_ext_dir = [-ones(Ng / 2, 1); ones(Ng / 2, 1)]; % applied constant torque on all particles
%     else
%         T_ext_dir = [-ones((Ng - 1) / 2, 1); 0; ones((Ng - 1) / 2, 1)];
%     end
%     theta_drive_init = rand(Ng, 1) * 2 * pi;
    theta_drive_init = (1.5 + (0:Ng-1)' / 5) * pi;
%     theta_drive_init = (1.5 + [0; 0.25; 0.75; 0.5; 0; 0.25; 0.75; 0.5]) * pi;
end
ext_idx1 = 2 * Ng + 1; % starting index for torque in the force&torque array
ext_idx2 = 3 * Ng; % ending index for torque in the force&torque array

xyt = xyt_init;
theta_drive = theta_drive_init;
t = 0;

nt_record = 100; % store particle positions every nt_record time step
xyt_record = zeros(3 * Ng, Nt / nt_record + 1); % array to store particle positions during MD simulation
xyt_record(:, 1) = xyt;

Vel = zeros(3 * Ng, 1); % translation and rotation velocity
Acc = Force_Bumpy(xyt, Epara); % force and torque calculation
Acc(ext_idx1:ext_idx2) = Acc(ext_idx1:ext_idx2) + T_drive * sin(f_drive * t + theta_drive); % add applied external torque
Acc_old = Acc;

% MD using Verlet method
for nt = 1:Nt
    % update position by a full step
    xyt = xyt + Vel * dt + Acc_old * dt_sq_half;

    if mod(nt, nt_record) == 0
        xyt_record(:, nt / nt_record + 1) = xyt;
    end
    if mod(nt, 100000) == 0
        fprintf('nt: %d\n', nt);
    end

    % update force by a full step
    t = t + dt;
    Acc = Force_Bumpy(xyt, Epara);
    Acc(ext_idx1:ext_idx2) = Acc(ext_idx1:ext_idx2) + T_drive * sin(f_drive * t + theta_drive);
    % correct force with damping term
    Acc = (Acc - B * (Vel + Acc_old * dt_half)) / B_denom;
    % update velocity by a full step
    Vel = Vel + (Acc + Acc_old) * dt_half;
    Acc_old = Acc;
end

%% Motion Test -- apply a constant velocity, without damping, for test only
% dt = 0.01; % time size for constant energy MD
% Nt = 1e5; % total number of time steps for MD
% dt_half = dt / 2;
% dt_sq_half = 0.5 * dt * dt;
% B = 0; % damping coefficient, set to 0 for constant energy MD
% B_denom = 1 + B * dt_half;
% 
% xyt = xyt_init;
% 
% nt_record = 10; % store particle positions every nt_record time step
% xyt_record = zeros(3 * Ng, Nt / nt_record + 1); % array to store particle positions during MD simulation
% xyt_record(:, 1) = xyt;
% 
% Vel = zeros(3 * Ng, 1); % translation and rotation velocity
% Vel(1:Ng) = 0.015 * [ones(Ng / 2, 1); -ones(Ng / 2, 1)];
% Acc = Force_Bumpy(xyt, Epara); % force and torque calculation
% Acc_old = Acc;
% 
% % MD using Verlet method
% for nt = 1:Nt
%     % update position by a full step
%     xyt = xyt + Vel * dt + Acc_old * dt_sq_half;
% 
%     if mod(nt, nt_record) == 0
%         xyt_record(:, nt / nt_record + 1) = xyt;
%     end
%     if mod(nt, 100000) == 0
%         fprintf('nt: %d\n', nt);
%     end
% 
%     % update force by a full step
%     Acc = Force_Bumpy(xyt, Epara);
%     % correct force with damping term
%     Acc = (Acc - B * (Vel + Acc_old * dt_half)) / B_denom;
%     % update velocity by a full step
%     Vel = Vel + (Acc + Acc_old) * dt_half;
%     Acc_old = Acc;
% end

%% make a movie of the MD simulation
dtheta = 2 * pi / Ns;

h = figure;
fig = gcf;
set(fig, 'DoubleBuffer', 'on');
set(gca, 'NextPlot', 'replace');

% opengl('software');
mf = 'Files/Roll';
mov = VideoWriter(mf, 'MPEG-4');
mov.FrameRate = 10;
open(mov);

for nt = 1:100:Nt/nt_record
    x = zeros(N, 1);
    y = zeros(N, 1);
    for ng = 1:Ng
        xc = xyt_record(ng, nt);
        yc = xyt_record(ng + Ng, nt);
        theta_s = xyt_record(ng + 2 * Ng, nt) + (0:Ns-1) * dtheta;
        idx1 = idx_start(ng);
        idx2 = idx_end(ng);
        x(idx1:idx2) = xc + Rg * cos(theta_s);
        y(idx1:idx2) = yc + Rg * sin(theta_s);
    end

    figure(h); hold on; box off;
    
    for ng = 1:Ng
        for i = idx_start(ng):idx_end(ng)
            rectangle('Position', [x(i) - Rb, y(i) - Rb, Db, Db],...
                      'Curvature', [1 1], 'EdgeColor', 'b', 'FaceColor', 'b');
        end
    end
    
    % plot([-1 2 * Ng], [0 0], 'k-', 'linewidth', 2)
    for ns = round(-10 / Db):round(2 * Ng / Db)
        rectangle('Position', [(ns - 1) * Db - Rb, -Rb, Db, Db],...
                  'Curvature', [1 1], 'EdgeColor', 'k', 'FaceColor', 'k');
    end
    
    axis equal
    axis off
    xlim([-10 16])
    ylim([-0.2 3.5])

    set(gcf, 'color', 'w'); 

    writeVideo(mov,getframe);
    clf;
end

close(h);
close(mov);
