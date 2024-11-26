% author: Mauro Morini
% last modified: 25.12.24
clc;clear;close all

% Intitializations
a = 0;
b = 2;
T = pi;
plotTimes = [0.1:0.02:T-0.1];
Hmax = 0.1;
Hmin = 0.01;
r = 1;
dt = Hmin*r;
projectionType = "";
meshTransformFrequency = 1;
meshCreationType = "refMid10";
meshTransformType = "addAndRemoveRand";

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


%% plot at time t
figure(1)
tPlot = 1;
[~, plotIdx] = min(abs(tPlot-plotTimes));
plotMesh = a:Hmin:b;
[Mesh, t, U] = PlotSol{plotIdx}.getSolution();
plot(plotMesh, uExact(plotMesh, t), 'Color',"#A2142F")
hold on
plot(Mesh.p, U,'Color', "#0072BD")
plot(Mesh.p,zeros(size(Mesh.p, 1)),'.', 'Color', "#77AC30")
hold off
legend("exact sol", "numerical sol")
title("T = " + t)

%% plot all
figure(2)
for i = 1:length(plotTimes)
    [Mesh, t, U] = PlotSol{i}.getSolution();
    plot(plotMesh, uExact(plotMesh, t), 'Color',"#A2142F")
    hold on
    plot(Mesh.p, U,'Color', "#0072BD")
    plot(Mesh.p,zeros(size(Mesh.p, 1)),'.', 'Color', "#77AC30")
    hold off
    ylim([-1,1])
    legend("exact sol", "numerical sol")
    title("T = " + t)
    drawnow
    pause(0.01)
end