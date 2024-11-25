% author: Mauro Morini
% last modified: 19.11.24
clc;clear;close all

% Intitializations
a = 0;
b = 2;
T = pi;
projectionType = "";
meshTransformFrequency = 1;
meshCreationType = "";
meshTransformType = "shiftHh";
Hmax = 2.^(-(2:0.5:4.5));
errors = zeros(1, length(Hmax));

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

for i = 1:length(Hmax)
    h = Hmax(i);
    dt = h/100;
    Mesh = createMeshRoutines([a,b],[h, h/100],meshCreationType);
    [UT, MeshT, Tend] = rothe1D(Mesh, v0, v1, f, dt, T, projectionType, meshTransformFrequency, meshTransformType);
    errors(1,i) = FEM1D.errorsLinear1D(MeshT.t,MeshT.p,UT,@(x)0,@(x) uExact(x,Tend));
end

% plots
figure(1)
tld = tiledlayout("flow");
nexttile
plot(MeshT.p, UT, Mesh.p, uExact(Mesh.p,T))
xlabel("x")
ylabel("y")
legend("u_{h}", "u_{exact}")
nexttile
loglog(Hmax,errors,Hmax, Hmax.^2, '--')
xlabel("Hmax")
ylabel("error")
ylim([0,1])
legend("L2err", "Hmax^2")