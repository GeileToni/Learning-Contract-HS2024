% author: Mauro Morini
% last modified: 25.11.24
function [EndTimeSol, PlotSol] = rothe1D(Mesh, v0, v1, f, dt, T, projectionType, meshTransformer, plotTimes)
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
% meshTransformer: MeshTransformer type containing information for the mesh
%                   transformation
% plotTimes: (1, N2) time vector with times in [0, T] at which the
%           solutions, the actual time and the Mesh should be output for
%           plots
%
% Outputs: 
% U: (NT,1) coefficient vecor for U at time T
% MT: Mesh1D object for the final mesh corresponding to UT
% T: scalar exact end time at which UT is numerically approximated, may be
%       different than original T

% Initializations
[p,~,t] = Mesh.getPet();
c = @(x) 1;
U = {v0(p)};
MeshList = {Mesh, Mesh};
time = 0:dt:T;

% adapt times to plot
if ~exist('plotTimes', 'var')
    plotTimes = [Inf];
else
    plotTimes = unique(plotTimes, 'sorted');
    plotTimes = plotTimes((dt*2 <= plotTimes) & (plotTimes <= T));
end
PlotSol = cell(1, length(plotTimes));
plotSolIdx = 1;

% get matrices for initial mesh
A = FEM1D.stiffnessMatrix1D(p,t,c);
F = FEM1D.loadVector1D(p, t, @(x) f(x,0));
M = FEM1D.massMatrix1D(p, t, c);
M = diag(sum(M, 2));                                % mass-lumped
intIdx = 2:(size(A,2)-1);

% Calculate U1
RHS = (dt^2*F(intIdx) - (A(intIdx,intIdx)*dt^2 - 2*M(intIdx,intIdx))*U{1}(intIdx) + 2*dt*M(intIdx,intIdx)*v1(p(intIdx)))/2;
U{2} = [0;M(intIdx,intIdx)\RHS;0];

% Leap Frog iteration
for i = 2:length(time)-1

    % change mesh
    [meshTransformer, Mesh, meshWasChanged] = meshTransformer.isTimeToChange(Mesh, time, i);
    MeshList{i+1} = Mesh;
    [p,~,t] = Mesh.getPet();
    
    if meshWasChanged       
        % Assemble matrices
        M = FEM1D.massMatrix1D(p, t, c);
        M = diag(sum(M, 2));                                % mass-lumped
        A = FEM1D.stiffnessMatrix1D(p,t,c);
    end
    
    % project previous U onto current mesh
    uPrev = project(U{i-1}, MeshList{i-1}, Mesh, projectionType);           
    uNow = project(U{i}, MeshList{i}, Mesh, projectionType);
    % uPrev = U{i-1};
    % uNow = U{i};

    F = FEM1D.loadVector1D(p, t, @(x) f(x,time(i+1)));
    intIdx = 2:(size(A,2)-1);

    % Solve system
    RHS = dt^2*F(intIdx) + (2*M(intIdx,intIdx) - dt^2*A(intIdx,intIdx))*uNow(intIdx) - M(intIdx,intIdx)*uPrev(intIdx);
    U{i+1} = [0;M(intIdx,intIdx)\RHS;0];

    % save Solution
    if plotSolIdx <= length(plotTimes) && abs(plotTimes(plotSolIdx) - time(i+1)) < dt/2
        WObj = WaveSolution(Mesh, time(i+1), U{i+1});
        PlotSol{plotSolIdx} = WObj;
        plotSolIdx = plotSolIdx + 1;
    end
end

% UT = U{end};
% MeshT = Mesh;
% T = t(end);

EndTimeSol = WaveSolution(Mesh, time(end), U{end});
end