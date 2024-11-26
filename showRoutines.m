% author: Mauro Morini
% last modified: 18.11.24
clc;clear;close all;

% Intitializations
Hmax = 0.1;
Hmin = 0.01;
a = 0;
b = 2;
T = pi;
meshCreationType = "refH/2";
meshTransformType = "addAndRemoveRand";
paramSpec = 4;

% create Mesh
Mesh = createMeshRoutines([a,b],[Hmax, Hmin],meshCreationType);
Mesh2 = changeMeshRoutines(Mesh, meshTransformType, paramSpec);

f1 = figure(1);
f2 = figure(2);
Mesh.plotMesh(f1);
Mesh2.plotMesh(f2);