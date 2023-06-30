function xyt_equ = Damped_MD_Minimization(xyt, Epara, Fthresh, dt, Nt)
%% Initialization
dt_half = dt / 2;
dt_sq_half = 0.5 * dt * dt;
B = 0.1;
B_denom = 1 + B * dt_half;

v_num = size(xyt, 1);
Vel = zeros(v_num, 1);
Acc = Force_Bumpy(xyt, Epara);
Acc_old = Acc;

if max(abs(Acc)) < Fthresh
    xyt_equ = xyt;
    return;
end
%% Damped MD
for nt = 1:Nt
    % MD using Verlet method
    % update position by a full step
    xyt = xyt + Vel * dt + Acc_old * dt_sq_half;
    % update force by a full step
    Acc = Force_Bumpy(xyt, Epara);

    % exit when forces are smaller than threshold
    if max(abs(Acc)) < Fthresh
        xyt_equ = xyt;
        return;
    end
    
    % correct force with damping term
    Acc = (Acc - B * (Vel + Acc_old * dt_half)) / B_denom;
    % update velocity by a full step
    Vel = Vel + (Acc + Acc_old) * dt_half;
    Acc_old = Acc;
end
%%
xyt_equ = xyt;