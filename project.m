% author: Mauro Morini
% last modified: 03.11.24
function u = project(uOld, MeshOld, MeshNew, projParam)
% projects solution uOld onto new mesh p
%
% Inputs:
% uOld: (N,1) vector with function coefficients as entries
% pOld: (N,1) pointvector of old mesh
% p: (N,1) pointvector of new mesh
% projParam: type of projection, default is linear
%
% Outputs:
% u: (N,1) coefficient vector projected onto mesh p

[pOld,~,t1] = MeshOld.getPet();
[p,~,t2] = MeshNew.getPet();

switch projParam
    case "L2"
        M = FEM1D.projMassMatrix1D(pOld,t1,p,t2);
        Mp = FEM1D.massMatrix1D(p',t2,@(x) 1);
        u = Mp\(M*uOld);
    otherwise 
        u = interp1(pOld, uOld, p);
end
end