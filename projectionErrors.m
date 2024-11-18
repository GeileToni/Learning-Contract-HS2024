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
L2err = zeros(2,length(H));
f = @(x) exp(x.^2);
seed = 3;

% refine one element and take out one point 
for i = 1:length(H)
    Mrng = Mesh1D([0, pi], [H(i), H(i)/100]);
    Mrng = Mrng.refine(H(i)/3, 2:(length(Mrng.p)/2));
    Mrng2 = Mrng.createRngMesh(seed);
    uOld = f(Mrng.p);
    uNew = project(uOld,Mrng.p,Mrng2.p,"L2");
    L2err(1,i) = FEM1D.L2ProjectionErrorLinear(uOld, Mrng.p, uNew, Mrng2.p);
end


figure
loglog(H, L2err(1,:), H, L2err(2,:), H, H.^2,'--')
legend("L2err rand","Hmax^2")
xlabel("Hmax")

% figure
% plot(Mrng.p, uOld, Mrng2.p,uNew)
% legend("uOld", "uNew")
% Mrng.plotMesh(figure)