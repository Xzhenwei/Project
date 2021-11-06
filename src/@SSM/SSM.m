classdef SSM < Manifold
    %SSM Spectral Submanifold Computation
    %   Detailed explanation goes here
    
    properties
       contOptions = cocoOptions();
       FRCOptions = FRCOptions();
       PSDOptions = PSDOptions();
       
       ssmSEulerTimeDisp = true;
       tol = 1e-9;
       w_l = [];
       X_l = [];
    end
    
    methods
        %% other methods
        
        FRC = extract_FRC(obj, parName, parRange, order)
                
        BB = extract_backbone(obj,modes,order,outdof)
        
        varargout = FRC_cont_ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        
        [FRC] = FRC_level_set(obj, resMode, order, parName, parRange)        
        
        
        [w,Gzz] = compute_ssmPSD(obj, PSDpair, W0, R0, method)
        
        [wss,Gss] = extract_PSD(obj, PSDpair, ORDER, method, freq_range, clusterRun)
        
        [w_l, X_l] = compute_analyticSSMPSD(obj,PSDpair,freq_range ,clusterRun)
    end
end

