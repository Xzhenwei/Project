classdef SSM < Manifold
    %SSM Spectral Submanifold Computation
    %   Detailed explanation goes here
    
    properties
       contOptions = cocoOptions();
       FRCOptions = FRCOptions();
    end
    
    methods
        %% other methods
        
        FRC = extract_FRC(obj, parName, parRange, order)
                
        BB = extract_backbone(obj,modes,order,outdof)
        
        varargout = FRC_cont_ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        
        [FRC] = FRC_level_set(obj, resMode, order, parName, parRange)        
        
        p = indirect_Euler_SSM(obj, N,T0,PSD,f,m,Wnode,R0)
        
        [w,Gzz] = compute_ssmPSD(obj, PSDpair, ORDER, method,clusterRun)
    end
end

