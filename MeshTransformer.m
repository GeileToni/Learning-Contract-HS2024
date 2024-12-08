% author: Mauro Morini
% last modified: 03.12.24
classdef MeshTransformer
    % contains information about the mesh change types, the frequency and
    % records how many changes have been made

    properties
        frequencyType           % String describing if the mesh change occurs every set step, at given indices or at given times
        frequency               % depending on frequencyType is either a single int, a index vector or a time vector
        meshTransformType       % String describing the type of mesh change, could possibly be a vector 
        transformCounter=0;     % counts how many times mesh changes have been implemented
    end

    methods
        function obj = MeshTransformer(frequencyType, frequency, meshTransformType)
            % standard constructor
            obj.frequencyType = frequencyType;
            obj.frequency = frequency;
            obj.meshTransformType = meshTransformType;
        end

        function [obj, Mesh] = changeMeshWithRoutine(obj, Mesh)
            % changes Mesh based on properties 
            meshTransformTypeIdx = mod(obj.transformCounter+1, length(obj.meshTransformType)) + 1;
            type = obj.meshTransformType(meshTransformTypeIdx);
            paramSpec = obj.transformCounter;   % parameter for change
            n = 3;                    % number of points removed/refined
            switch type
                case "rng"
                    Mesh = Mesh.createRngMesh(paramSpec);
                case "shiftH/4"
                    Mesh = Mesh.shiftMesh(Mesh.Hmax/4);
                case "shiftHh"
                    Mesh = Mesh.shiftMesh(Mesh.Hmax/4+Mesh.Hmin);
                case "shiftBackAndForth"
                    hShift = (Mesh.Hmax/4+Mesh.Hmin)*(-1)^(paramSpec);
                    Mesh = Mesh.shiftMesh(hShift);
                case "removeRand1"
                    % removes 1 random points
                    Mesh = Mesh.removeRand(paramSpec);
                case "removeRand"
                    % removes n random points
                    Mesh = Mesh.removeRand(n, paramSpec);
                case "addAndRemoveRand"
                    rng(paramSpec);
                    for i = 1:n
                        refinableTIdx = Mesh.getRefinableElements(2);
                        rdmIdx = randi([1,length(refinableTIdx)],1,1);
                        Mesh = Mesh.refineByFact(refinableTIdx(rdmIdx));
                        Mesh = Mesh.removeRand(rdmIdx);
                    end
                otherwise % no change
            end
        end

        function [obj, Mesh, meshWasChanged] = isTimeToChange(obj, Mesh, timeVec, currentIdx)
            % depending on frequency checks if this iteration/ at this time
            % there should be made a mesh change and treats it accordingly
            meshWasChanged = false;
            switch obj.frequencyType
                case "regular"
                    if ~(mod(currentIdx-1, obj.frequency) == 0)
                        return
                    end
                case "time"
                    dt = abs(timeVec(1)-timeVec(2));
                    itsTime = abs(obj.frequency - timeVec(currentIdx)) < dt/2;
                    if ~any(itsTime)
                        return
                    end
                otherwise % do nothing
                    return
            end
            [obj, Mesh] = obj.changeMeshWithRoutine(Mesh);
            obj.transformCounter = obj.transformCounter + 1;
            meshWasChanged = true;
        end
    end
end