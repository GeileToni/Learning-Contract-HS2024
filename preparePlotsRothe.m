% author: Mauro Morini
% last modified: 09.12.24
clc;clear;close all

MeshCreationTypes = ["rng1", "rngRef5,1", "refH/2", "refMid10", "refInnerAllSec10", "refInnerAll10",""];
MeshTransformTypes = ["rng", "shiftH/4", "shiftHh", "shiftBackAndForth", "removeRand1", "removeRand", "addAndRemoveRand"];

% save rates 
saveIdx = [1, 3, 4, 7];
solutions = cell(2, 8);
for k1 = 1:2
for k3 = 1:2
for k2 = 1:length(saveIdx)

% Intitializations
a = 0;
b = 2;
T = pi;
projectionType = "L2";
meshTransformFrequency = k3;
frequencyType = "regular";
meshCreationType = MeshCreationTypes(3);
meshTransformType = MeshTransformTypes(saveIdx(k2));
MT = MeshTransformer(frequencyType, meshTransformFrequency, meshTransformType);
Hmax = 2.^(-(2:0.5:6));
errors = zeros(2, length(Hmax));
r = 0.9;
Dof = k1;

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

% Calculation of errors

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

solutions{k1,(k3-1)*4 + k2} = struct('errors', errors, 'MeshTransformer', MT, 'Dof', Dof, 'Hmax', Hmax, 'EndSol', EndSol);

end
end
end
save("errorsAndSolData","solutions")


%% PLOT
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

for k2 = 1:size(solutions,2)
% Extract solutions
Dof = 2;
idx = k2;
errors = solutions{Dof, idx}.errors;
[MeshT, Tend, UT] = solutions{Dof, idx}.EndSol.getSolution();
Hmax = solutions{Dof, idx}.Hmax;
meshTransformType = solutions{Dof, idx}.MeshTransformer.meshTransformType;
frequency = solutions{Dof, idx}.MeshTransformer.frequency;

% plots
p = MeshT.getPet();
figure()
tld = tiledlayout("flow");
nexttile
plot(p, uExact(p,Tend),p, UT)
xlabel("x")
ylabel("y")
legend("u_{exact}","u_{h}")
nexttile
loglog(Hmax,errors(1,:),Hmax, Hmax.^2, '--',Hmax, Hmax, '--', Hmax, errors(2,:), Hmax, Hmax.^3,'--', Hmax, Hmax.^4,'--')
xlabel("Hmax")
ylabel("error")
% ylim([0,1])
legend("L2err", "Hmax^2", "Hmax", "H1err", "Hmax^3","Hmax^4")
title(tld, "Dof: " + Dof + "| mesh Transform type = " + meshTransformType + "| frequency = " + frequency)

end

