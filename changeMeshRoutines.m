% author: Mauro Morini
% last modified: 19.11.24
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
        case "shiftBackAndForth"
            hShift = (M.Hmax/4+M.Hmin)*(-1)^(paramSpec);
            M = M.shiftMesh(hShift);
        case "removeRand1"
            % removes 1 random points
            M = M.removeRand(paramSpec);
        case "removeRand"
            % removes n random points
            n = 3;
            M = M.removeRand(n, paramSpec);
        case "addAndRemoveRand"
            n = 5;
            rng(paramSpec);
            for i = 1:n
                refinableTIdx = M.getRefinableElements(2);
                rdmIdx = randi([1,length(refinableTIdx)],1,1);
                M = M.refineByFact(refinableTIdx(rdmIdx));
                M = M.removeRand(rdmIdx);
            end
        otherwise % no change
    end
end