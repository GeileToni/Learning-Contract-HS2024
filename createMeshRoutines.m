% author: Mauro Morini
% last modified: 03.11.24
function M = createMeshRoutines(ab, H, type)
% Contains mesh creation routines 
M = Mesh1D(ab, H);
switch type(1)
    case "rng1"
        % create rng mesh
        M = M.createRngMesh(1);
    case "rngRef5,1"
        % refine one random element by factor 5 with seed 1
        r = rng(1, "twister");
        nT = size(M.t,1);
        numRef = rand(1);
        numRef = floor(numRef*(nT-1));
        refIdx = rand(numRef,1);
        refIdx = floor(refIdx*(nT-1)+1);
        refIdx = unique(refIdx, "sorted");
        Href = M.Hmax/5;
        M = M.refine(Href, refIdx);
    case "refH/2"
        M = M.refine(M.Hmax/2, 1:size(M.t,1));
    case "refMid10"
        m = (M.a + M.b)/2;
        containsm = find(M.p(M.t(:,1)) <= m & m <= M.p(M.t(:,2)));
        M = M.refine(M.Hmax/10, containsm);
    case "refInnerAllSec10"
        refIdx = 2:2:size(M.t,1)-1;
        M = M.refine(M.Hmax/10, refIdx);
    case "refInnerAll10"
        refIdx = 2:size(M.t,1)-1;
        M = M.refine(M.Hmax/10, refIdx);
    otherwise
end