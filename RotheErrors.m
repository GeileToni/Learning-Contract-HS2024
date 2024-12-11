% author: Mauro Morini
% last modified: 03.12.24
clc;clear;close all

MeshCreationTypes = ["rng1", "rngRef5,1", "refH/2", "refMid10", "refInnerAllSec10", "refInnerAll10",""];
MeshTransformTypes = ["rng", "shiftH/4", "shiftHh", "shiftBackAndForth", "removeRand1", "removeRand", "addAndRemoveRand"];

% Intitializations
a = 0;
b = 2;
T = pi;
projectionType = "L2";
meshTransformFrequency = 2;
frequencyType = "regular";
meshCreationType = MeshCreationTypes(3);
meshTransformType = MeshTransformTypes(1);
MT = MeshTransformer(frequencyType, meshTransformFrequency, meshTransformType);
Hmax = 2.^(-(2:5));
errors = zeros(2, length(Hmax));
r = 0.9;
Dof = 2;

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
    dt = r*h/10/(Dof^2);
    Mesh = createMeshRoutines([a,b],[h, h/10],meshCreationType);
    Mesh.r = Dof;
    EndSol = rothe1D(Mesh, v0, v1, f, dt, T, projectionType, MT);
    [MeshT, Tend, UT] = EndSol.getSolution();
    [p,~,t] = MeshT.getPet();
    [errors(1,i), errors(2,i)] = FEM1D.errors1D(t,p,UT,@(x) dudx(x,Tend),@(x) uExact(x,Tend));
end

%% plots
p = MeshT.getPet();
figure(1)
tld = tiledlayout("flow");
nexttile
plot(p, uExact(p,T),p, UT)
xlabel("x")
ylabel("y")
legend("u_{exact}","u_{h}")
nexttile
loglog(Hmax,errors(1,:),Hmax, Hmax.^2, '--',Hmax, Hmax, '--', Hmax, errors(2,:), Hmax, Hmax.^3,'--', Hmax, Hmax.^4,'--')
xlabel("Hmax")
ylabel("error")
% ylim([0,1])
legend("L2err", "Hmax^2", "Hmax", "H1err", "Hmax^3","Hmax^4")
title(tld, "Dof: " + Dof + "| Projection type: " + projectionType + "| r = " + r + "| mesh Transform type = " + meshTransformType + "| frequency = " + meshTransformFrequency)


%% Calculate convergence
log_h = log(Hmax(3:end));
R = log_h'.^(0:1);
% convergence order L2 error
log_L2_err = log(errors(1,3:end));
coeff_L2_err = R\log_L2_err';
fprintf('convergence order L2 error: %1.4f\n',coeff_L2_err(2));
% convergence order H1 error 
log_H1_err = log(errors(2,3:end));
coeff_H1_err = R\log_H1_err';
fprintf('convergence order H1 error: %1.4f\n',coeff_H1_err(2));