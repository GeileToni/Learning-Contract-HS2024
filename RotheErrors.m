% author: Mauro Morini
% last modified: 03.12.24
clc;clear;close all

MeshCreationTypes = ["rng1", "rngRef5,1", "refH/2", "refMid10", "refInnerAllSec10", "refInnerAll10",""];
MeshTransformTypes = ["rng", "shiftH/4", "shiftHh", "shiftBackAndForth", "removeRand1", "removeRand", "addAndRemoveRand",""];

% Intitializations
a = 0;
b = 2;
T = pi;
projectionType = "";
meshTransformFrequency = 1;
frequencyType = "regular";
meshCreationType = MeshCreationTypes(3);
meshTransformType = MeshTransformTypes(4);
Hmax = 2.^(-(2:0.5:6));
errors = zeros(2, length(Hmax));
r = 0.9;

% functions
syms x t
uExact = sin(pi*x)*cos(pi*t);
v1 = diff(uExact, t, 1);
dudx = diff(uExact, x, 1);
f = diff(uExact, t, 2) - diff(uExact, x, 2);
uExact = matlabFunction(uExact, "Vars", [x, t]);
dudx = matlabFunction(dudx, "Vars", [x, t]);
v1 = matlabFunction(v1, "Vars", [x, t]);
f = matlabFunction(f, "Vars", [x, t]);
v0 = @(x) uExact(x, 0);
v1 = @(x) v1(x,0);

for i = 1:length(Hmax)
    h = Hmax(i);
    dt = r*h/10;
    Mesh = createMeshRoutines([a,b],[h, h/10],meshCreationType);
    MT = MeshTransformer(frequencyType, meshTransformFrequency, meshTransformType);
    EndSol = rothe1D(Mesh, v0, v1, f, dt, T, projectionType, MT);
    [MeshT, Tend, UT] = EndSol.getSolution();
    [errors(1,i), errors(2,i)] = FEM1D.errorsLinear1D(MeshT.t,MeshT.p,UT,@(x) dudx(x,Tend),@(x) uExact(x,Tend));
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
loglog(Hmax,errors(1,:),Hmax, Hmax.^2, '--',Hmax, Hmax, '--', Hmax, errors(2,:))
xlabel("Hmax")
ylabel("error")
% ylim([0,1])
legend("L2err", "Hmax^2", "Hmax", "H1err")