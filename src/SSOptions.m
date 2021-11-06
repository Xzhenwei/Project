classdef SSOptions < matlab.mixin.SetGet
    %SSOptions Options for DynamicalSystems class
    properties
        ssMethod
        outDOF = []; % output degree of freedom
        tol
        
    end
    methods
        
        function set.ssMethod(obj,ssMethod)
            switch lower(ssMethod)
                case 'direct'
                    obj.ssMethod = 'direct';
                case 'indirect'
                    obj.ssMethod = 'indirect';
                otherwise
                    error('Unknown stochastic input type: set direct or indirect method types')
            end
        end
        
        function set.tol(obj,Tol)
            obj.tol = Tol;
        end
        
    end
    
end