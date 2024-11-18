% author: Mauro Morini
% last modified: 09.11.24
clc; clear; close all;
% Testscript 

%% refine test
M = Mesh1D([0, 1], [0.1, 0.001]);
M = M.refine(M.Hmax/2, 1:size(M.t,1));
% plot
f = figure(1);
M.plotMesh(f);

%% shift test
M = Mesh1D([-2,2], [0.4, 0.001]);
M = M.shiftMesh(0.06);
f = figure(2);
M.plotMesh(f);

%% coarsen test
M = Mesh1D([0, 1], [0.1, 0.001]);
M = M.refine(M.Hmax/10, [1]);
for i = 1:8
    M = M.removePoints([2]);
end
f = figure(3);
M.plotMesh(f);

%% project

M1 = Mesh1D([0, pi], [0.3, 0.001]);
x = 1;
M1 = M1.refine(M1.Hmax/2, M1.findElementAt(x));
M2 = M1.removePoints(M1.findPointClosestTo(x));
f = figure(3);
M2.plotMesh(f);
uOld = sin(M1.p);
uNew = project(uOld,M1.p,M2.p,"L2");
hold on
plot(M1.p, uOld, '-o', M2.p, uNew, '-x')
hold off

%% projection mass matrix 
M = Mesh1D();
M1 = FEM1D.massMatrix1D((M.p)',M.t,@(x) 1);
M2 = FEM1D.projMassMatrix1D(M.p,M.t,M.p,M.t);
max(abs(M1-M2),[],"all")

%% Test random mesh
seed = 2;
M = Mesh1D();
M = M.createRngMesh(seed);
f = figure(4);
M.plotMesh(f)

%% Test create Mesh routine rng
a = 0;
b = 1;
Hmax = 0.1;
Hmin = Hmax/100;
M = createMeshRoutines([a,b],[Hmax, Hmin], "rngRef5,1");
f = figure(4);
M.plotMesh(f);

%% Test shift
a = 0;
b = 2;
M = createMeshRoutines([a,b], [0.1, 0.001],"");
M = changeMeshRoutines(M, "shiftHh");

%% Test create Mesh routine refMid10
a = 0;
b = 1;
Hmax = 0.1;
Hmin = Hmax/100;
M = createMeshRoutines([a,b],[Hmax, Hmin], "refMid10");
f = figure(4);
M.plotMesh(f);
