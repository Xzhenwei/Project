classdef StochasticSystem < DynamicalSystem
    
    
    properties

        filterPSD = []      % input PSD defined by Linear SYSTEM
        samplePSD = []      % input PSD defined by sampled PSD
        input
        inputForcingType
        linear              % if model is linear
        nRealization 
        timeSpan
        nPoints
        forcingdof % forcing Degree Of Freedom
        
        timevector = []
        
        gFactor = 1; % this is the factor used in indirect process generation
        Fsto = [] % an vector of stochastic forcing realization
        
        sdeImpTimeDisp = false;
        
        SSOptions = SSOptions()
        
    end
    


    methods

        
        %% other methods
        function add_random_forcing(obj,nRealization,timeSpan,nPoints,forcingdof)
            obj.nRealization = nRealization;
            obj.timeSpan = timeSpan;
            obj.nPoints = nPoints;
            obj.forcingdof = forcingdof; % forcing DOF
            obj.timevector = [];
        end
        
        
        [r, drdqdd,drdqd,drdq, c0] = residual(obj, q, qd, qdd, t)
        
        %%% stochastic methods
        [F_sto] = generate_stochastic(obj)
        fst = compute_fstochastic(obj,t)
        
        [X,V]=forward_Heun(obj,X0,V0,N,T0,PSD)
        
	    [X,V]=implicit_Mid_Point(obj,N,T0,PSD)

        [w,PHI]=sde_solver(obj,SDEmethod,PowerSpectralPair)
        
        [w,linear_analytic]=compute_linear_PSD(obj,PSDpair,freq_range)
        
        [w,PSD] = monte_carlo_average(obj,method,PSDpair,nRealization)
        
        [w,Gz] = galerkin_proj (obj, V, PSD, SDEmethod, PSDpair)
        
        [w_galerkin, PSD_galerkin] = monte_carlo_galerkin(obj, method, PSDpair)
        
%         input_PSD(obj)

    end
end

