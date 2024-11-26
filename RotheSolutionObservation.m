% author: Mauro Morini
% last modified: 25.12.24
clc;clear;close all

% Intitializations
a = 0;
b = 2;
T = pi;
plotTimes = [0.1:0.2:T-0.1];
Hmax = 0.1;
Hmin = 0.01;
r = 0.4;
dt = Hmin*r;
projectionType = "";
meshTransformFrequency = 10;
meshCreationType = "refH/2";
meshTransformType = "shiftBackAndForth";

% functions
syms x t
uExact = sin(pi*x)*cos(pi*t);
v1 = diff(uExact, t, 1);
f = diff(uExact, t, 2) - diff(uExact, x, 2);
uExact = matlabFunction(uExact, "Vars", [x, t]);
v1 = matlabFunction(v1, "Vars", [x, t]);
f = matlabFunction(f, "Vars", [x, t]);
v0 = @(x) uExact(x, 0);
v1 = @(x) v1(x,0);

% calculate solution
Mesh = createMeshRoutines([a,b],[Hmax, Hmin],meshCreationType);
[EndSol, PlotSol] = rothe1D(Mesh, v0, v1, f, dt, T, projectionType, meshTransformFrequency, meshTransformType, plotTimes);


%% plot
plotIdx = 15;
plotMesh = a:Hmin:b;
[Mesh, t, U] = PlotSol{plotIdx}.getSolution();
plot(plotMesh, uExact(plotMesh, t), Mesh.p, U, Mesh.p, zeros(size(Mesh.p, 1)), '.')
legend("exact sol", "numerical sol")
title("T = " + t)

