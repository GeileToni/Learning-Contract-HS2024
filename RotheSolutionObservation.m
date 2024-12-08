% author: Mauro Morini
% last modified: 25.12.24
clc;clear;close all

MeshCreationTypes = ["rng1", "rngRef5,1", "refH/2", "refMid10", "refInnerAllSec10", "refInnerAll10",""];
MeshTransformTypes = ["rng", "shiftH/4", "shiftHh", "shiftBackAndForth", "removeRand1", "removeRand", "addAndRemoveRand"];

% Intitializations
a = 0;
b = 2;
T = pi;
Hmax = 2^(-4.5);
Hmin = Hmax/10;
r = 0.5;
Dof = 2;
r = r/Dof;
dt = Hmin*r;

% % possibly calculate stable dt
% pTemp = (a:Hmin:b)';
% tTemp = [(1:length(pTemp)-1)', (2:length(pTemp))'];
% ATemp = FEM1D.stiffnessMatrix1D(pTemp, tTemp, @(x) 1);
% MTemp = FEM1D.massMatrix1D(pTemp, tTemp, @(x) 1);
% dtStable = sqrt(eigs(MTemp,1,0)/eigs(ATemp,1));
% dt = dtStable;

plotTimes = [0.1:0.05:T-0.1];
%plotTimes = 2*dt:dt:T;
projectionType = "L2";
meshTransformFrequency = 2;
frequencyType = "regular";
meshCreationType = MeshCreationTypes(3);
meshTransformType = MeshTransformTypes(end);

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
MT = MeshTransformer(frequencyType, meshTransformFrequency, meshTransformType);
Mesh = createMeshRoutines([a,b],[Hmax, Hmin],meshCreationType);
Mesh.r = Dof;
[EndSol, PlotSol] = rothe1D(Mesh, v0, v1, f, dt, T, projectionType, MT, plotTimes);


%% plot at time t
figure(1)
tPlot = 3;
[~, plotIdx] = min(abs(tPlot-plotTimes));
plotMesh = a:Hmin:b;
[Mesh, t, U] = PlotSol{plotIdx}.getSolution();
p = Mesh.getPet();
plot(plotMesh, uExact(plotMesh, t), 'Color',"#A2142F")
hold on
plot(p, U,'Color', "#0072BD")
plot(p,zeros(size(p, 1)),'.', 'Color', "#77AC30")
hold off
legend("exact sol", "numerical sol")
title("T = " + t)

%% plot all
figure(2)
plotMesh = a:Hmin:b;
for i = 1:length(plotTimes)
    [Mesh, t, U] = PlotSol{i}.getSolution();
    p = Mesh.getPet();
    plot(plotMesh, uExact(plotMesh, t), 'Color',"#A2142F")
    hold on
    plot(p, U,'Color', "#0072BD")
    plot(p,zeros(size(p, 1)),'.', 'Color', "#77AC30")
    hold off
    ylim([-1,1])
    legend("exact sol", "numerical sol")
    title("T = " + t)
    drawnow
    pause(0.01)
end