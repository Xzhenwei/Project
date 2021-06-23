function fst = compute_fstochastic(obj,t)
% COMPUTE_FEXT We compute the external force at a particular time t 
% in a second-order mechanical system. 

T = obj.timevector;
Fext = obj.Fsto;

fst = interp1(T,Fext',t)';
end
