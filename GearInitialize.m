function xyt = GearInitialize(pid, Epara)
%%
rng(pid);

Ng = Epara.Ng;
Ns = Epara.Ns;
Rg = Epara.Rg;
Db = Epara.Db;

xyt = zeros(3 * Ng, 1);

xyt(1:Ng) = (0:Ng-1)' * 4 * Rg;
xyt(Ng+1:2*Ng) = Rg + Db;
xyt(2*Ng+1:3*Ng) = 2 * pi / Ns * rand(Ng, 1);