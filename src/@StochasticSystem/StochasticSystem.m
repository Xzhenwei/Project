classdef StochasticSystem < DynamicalSystem
    
    
    properties
%         System      % dynamical system object
%         SSM         % SSM object

        filterPSD = []      % input PSD defined by Linear SYSTEM
        samplePSD = []      % input PSD defined by sampled PSD
        input
        
        linear              % if model is linear
        nRealization 
        timeSpan
        nPoints
        forcingdof % forcing Degree Of Freedom
        
        timevector = []
        
        Fsto = [] % an vector of stochastic forcing realization
        
        sdeImpTimeDisp = true;
        
        SSOptions = SSOptions()
    end
    


    methods
        %% Constructor
%         function obj = StochasticSystem(ClassInput,ClassType)
%             %STOCHASTIC Construct an instance of this class
%             %   Detailed explanation goes here
%             switch lower(ClassType)
%                 case 'system'
%                     obj.System = ClassInput;        
%                 case 'ssm'
%                     obj.SSM = ClassInput;
%             end
%         end

        %% SET methods

        
%         %% GET methods        
%         function N = get.dimSystem(obj)
%             N = obj.System.N;
%         end
%         

        
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
        
        linear_analytic=compute_linear_PSD(obj,omega,PSD)
        
        [w,PSD] = monte_carlo_average(obj,method,PSDpair,nRealization,clusterRun)
        
        input_PSD(obj)
        %%% SSM
%         PSD = extract_PSD(obj, parRange, order, method)
%         
%         p=indirect_Euler_SSM(obj, nPoints, timeSpan, PSD, dimFilter, dimManifold, Wnode, R0)
    end
end

