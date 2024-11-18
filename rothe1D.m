% author: Mauro Morini
% last modified: 09.11.24
function [UT, MeshT, T] = rothe1D(Mesh, v0, v1, f, dt, T, projectionType, meshTransformFrequency, meshTransformType)
% Uses the rothe method for solving the 1d wave equation u_xx = u_tt with
% zero dirichlet b.c.
% Discretizes in time first then in space, calculates sequence fully discrete
% coefficients of solution u in a changing mesh, projects previous
% solutions into current space using a given projection method then
% calculates newest coefficients using Leap-Frog for given time step and
% given initial cond. with mass lumping.
% 
%
% Inputs:
% M: Mesh1D object containing initial mesh information
% v0: @(x) function handle for initial cond u(x,0) = v0;
% v1: @(x) function handle for initial cond du/dt(x,0) = v1;
% f: @(x,t) function handle for inhomogeneous part
% dt: scalar constant time step
% T: scalar end time
% projectionType: string defining the projection type
% meshTransformationFrequency: positive integer denoting denoting frequency
%           of mesh transformation (ex freq = 1 => every iteration mesh is changed
% meshTransformType: string denoting type of mesh transformation
%
% Outputs: 
% U: (NT,1) coefficient vecor for U at time T
% MT: Mesh1D object for the final mesh corresponding to UT
% T: scalar exact end time at which UT is numerically approximated, may be
%       different than original T

% Initializations
c = @(x) 1;
U = {v0(Mesh.p)};
MeshList = {Mesh, Mesh};
t = 0:dt:T;

% get matrices for initial mesh
A = FEM1D.stiffnessMatrix1D(Mesh.p,Mesh.t,c);
F = FEM1D.loadVector1D(Mesh.p, Mesh.t, @(x) f(x,0));
M = FEM1D.massMatrix1D(Mesh.p, Mesh.t, c);
M = diag(sum(M, 2));                                % mass-lumped
intIdx = 2:(size(A,2)-1);

% Calculate U1
RHS = (dt^2*F(intIdx) - (A(intIdx,intIdx)*dt^2 - 2*M(intIdx,intIdx))*U{1}(intIdx) + 2*dt*M(intIdx,intIdx)*v1(Mesh.p(intIdx)))/2;
U{2} = [0;M(intIdx,intIdx)\RHS;0];

% Leap Frog iteration
for i = 2:length(t)-1

    % change mesh
    if mod(i-1, meshTransformFrequency) == 0
        Mesh = changeMeshRoutines(Mesh, meshTransformType, i);
    end
    MeshList{i+1} = Mesh;

    % project previous U onto current mesh
    uPrev = project(U{i-1}, MeshList{i-1}.p, Mesh.p, projectionType);
    uNow = project(U{i}, MeshList{i}.p, Mesh.p, projectionType);
    
    % Assemble matrices
    M = FEM1D.massMatrix1D(Mesh.p, Mesh.t, c);
    M = diag(sum(M, 2));                                % mass-lumped
    A = FEM1D.stiffnessMatrix1D(Mesh.p,Mesh.t,c);
    F = FEM1D.loadVector1D(Mesh.p, Mesh.t, @(x) f(x,t(i+1)));
    intIdx = 2:(size(A,2)-1);

    % Solve system
    RHS = dt^2*F(intIdx) + (2*M(intIdx,intIdx) - dt^2*A(intIdx,intIdx))*uNow(intIdx) - M(intIdx,intIdx)*uPrev(intIdx);
    U{i+1} = [0;M(intIdx,intIdx)\RHS;0];
end

UT = U{end};
MeshT = Mesh;
T = t(end);

%u1 = uExact(Mesh.p,t(i+1));
%uExact = @(x,t)cos(t.*pi).*sin(x.*pi)
%plot(Mesh.p, u1, Mesh.p, U{i+1}, 'x')


end