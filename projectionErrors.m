% author: Mauro Morini
% last modified: 31.10.24
clc; clear; close all;

% Initializations
a = 0;
b = 1;
H = 2.^(-(1:7));
L2err = zeros(2,length(H));
f = @(x) exp(x.^2);

% refine one element and take out one point 
for i = 1:length(H)
    M1 = Mesh1D([0, pi], [H(i), H(i)/100]);
    x = 1;
    M1 = M1.refine(M1.Hmax/2, M1.findElementAt(x));
    M2 = M1.removePoints(M1.findPointClosestTo(x));
    uOld = f(M1.p);
    uNew = project(uOld,M1.p,M2.p,"L2");
    L2err(1,i) = FEM1D.L2ProjectionErrorLinear(uOld, M1.p, uNew, M2.p);
end

% refine one element and take out one point 
for i = 1:length(H)
    M1 = Mesh1D([0, pi], [H(i), H(i)/100]);
    M2 = M1.shiftMesh(3*H(i)/50);
    uOld = f(M1.p);
    uNew = project(uOld,M1.p,M2.p,"L2");
    L2err(2,i) = FEM1D.L2ProjectionErrorLinear(uOld, M1.p, uNew, M2.p);
end

figure
loglog(H, L2err(1,:), H, L2err(2,:), H, H.^2,'--')
legend("L2err coars", "L2err shift","Hmax^2")
xlabel("Hmax")

%%
clc; clear; close all;

% Initializations
a = 0;
b = 1;
H = 2.^(-(1:7));
projType = "";
meshCreationType = "";
meshChangeType = "rng";
L2err = zeros(2,length(H));
f = @(x) exp(x.^2);
seed = 3;

% refine one element and take out one point 
for i = 1:length(H)
    Mesh1 = createMeshRoutines([a,b], [H(i), H(i)/100], meshCreationType);
    Mesh2 = changeMeshRoutines(Mesh1,meshChangeType,3);
    uOld = f(Mesh1.p);
    uNew = project(uOld,Mesh1.p,Mesh2.p,projType);
    L2err(1,i) = FEM1D.L2ProjectionErrorLinear(uOld, Mesh1.p, uNew, Mesh2.p);
end


figure
loglog(H, L2err(1,:), H, L2err(2,:), H, H.^2,'--')
legend("L2err","Hmax^2")
xlabel("Hmax")

% figure
% plot(Mesh1.p, uOld, Mesh2.p,uNew,'--')
% legend("uOld", "uNew")
% Mesh1.plotMesh(figure)