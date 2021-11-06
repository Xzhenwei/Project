function [w,Gz] = GalerkinProj (V, SS, M, C, K)

Mt = V'*M*V; Ct = V'*C*V; Kt = V'*K*V;
tol = 1e-8;
SS.Fsto = SS.generate_stochastic();   n = SS.nPoints;
T0 = SS.timeSpan;
detT = T0/n; t_span = 0:detT:T0; 
        %     %% ODE45 solver
                f_p = @(t,q) Galerkin(SS,V,Mt,Ct,Kt,t,q);
                opts = odeset('RelTol',tol,'AbsTol',tol);
                [t_45,q_45] = ode45(f_p,t_span ,zeros(2,1), opts);
                q_45 = q_45';
                x = V.*q_45(1,:);

                [w,Gz] = crossPSDestimator(x(1,:),x(1,:),t_45);
                
    function dq = Galerkin(SS,V,Mt,Ct,Kt,t,q)
        dq1 = q(2) ;
        dq2 = (V'*SS.compute_fstochastic(t) - V'*SS.compute_fnl(V*q(1),V*q(2)) - Kt*q(1) - Ct*q(2))/Mt ;
        dq = [dq1;dq2];
    end
                
end


function [w,Gxy]=crossPSDestimator(x,y,t)
%%% Reference: Preunont Andre, Random vibration and Spectral Analysis Chp 12.9
%%%% numebr of sampling points 
N = length(x); 
%%%% time period of generation
T0 = max(t);
%%%% sampling frequency, time spacing
f0 = 1/T0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cx = fft(x)/N; Cy = fft(y)/N; Cys = conj(Cy);

Gxy = T0*(Cx.*Cys)/2/pi;
w = (1:length(Cx))*f0*2*pi;

end