% author: Mauro Morini
% last modified: 18.11.24
clc;clear;close all;

% Intitializations
Hmax = 0.5;
Hmin = Hmax/100;
a = 0;
b = 2;
T = pi;
meshTransformType = "rng";
meshCreationType = "";
paramSpec = 1;

% create Mesh
Mesh = createMeshRoutines([a,b],[Hmax, Hmin],meshCreationType);
Mesh2 = changeMeshRoutines(Mesh, meshTransformType, paramSpec);

f = figure(1);
hold on
Mesh.plotMesh(f);
Mesh2.plotMesh(f);
hold off