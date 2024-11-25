% author: Mauro Morini
% last modified: 25.11.24
classdef WaveSolution
    % small object containing information of a certain numerical solution at a
    % time T and it's current mesh 

    properties
        T               % Time at which numerical sol is recorded
        U               % (nP, 1) point array
        Mesh            % Mesh1D object 
    end

    methods
        function obj = WaveSolution(Mesh, T, U)
            obj.Mesh = Mesh;
            obj.T = T;
            obj.U = U;
        end

        function [Mesh, T, U] = getSolution(obj)
            Mesh = obj.Mesh;
            T = obj.T;
            U = obj.U;
        end
    end
end