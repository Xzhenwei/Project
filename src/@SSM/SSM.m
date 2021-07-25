classdef SSM < Manifold
    %SSM Spectral Submanifold Computation
    %   Detailed explanation goes here
    
    properties
       contOptions = cocoOptions();
       FRCOptions = FRCOptions();
       
       ssmSEulerTimeDisp = true;
    end
    
    methods
        %% other methods
        
        FRC = extract_FRC(obj, parName, parRange, order)
                
        BB = extract_backbone(obj,modes,order,outdof)
        
        varargout = FRC_cont_ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        
        [FRC] = FRC_level_set(obj, resMode, order, parName, parRange)        
        
        p = ssm_Euler_solver(obj, N,T0,PSD,f,m,Wnode,R0)
        
        [w,Gzz] = compute_ssmPSD(obj, PSDpair, W0, R0, method)
        
        [wss,Gss] = extract_PSD(obj, PSDpair, ORDER, method,clusterRun)
        
        [w, X_l] = compute_analyticPSD(obj,PSDpair)
    end
end

