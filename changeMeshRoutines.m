% author: Mauro Morini
% last modified: 09.11.24
function M = changeMeshRoutines(M, type, paramSpec)
% Changes a given input mesh based on a string declaring the type of change
% and an additional parameter for specification

switch type
    case "rng"
        M = M.createRngMesh(paramSpec);
    case "shiftH/4"
        M = M.shiftMesh(M.Hmax/4);
    case "shiftHh"
        M = M.shiftMesh(M.Hmax/4+M.Hmin);
    case "removeRand"
        % removes random point

    otherwise % no change
end
end