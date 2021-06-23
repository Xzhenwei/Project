function [ r, drdqdd,drdqd,drdq, c0] = residual_coupled(obj, q, qd, qdd, t )
%  RESIDUAL This function computes the residual needed for time integration 
% of
% 
% second-order system
% 
% $\mathbf{M}\ddot{\mathbf{q}} + \mathbf{C}\dot{\mathbf{q}} + \mathbf{F}(\mathbf{q}) 
% =\mathbf{F}_{ext}(t)$,
% 
% where we use the residual is defined as
% 
% $\mathbf{r}(\ddot{\mathbf{q}},\dot{\mathbf{q}},\mathbf{q}) = \mathbf{M}\ddot{\mathbf{q}} 
% + \mathbf{C}\dot{\mathbf{q}} + \mathbf{F}(\mathbf{q}) - \mathbf{F}_{ext}(t)$.

assert(obj.order == 2, ' residual can only be computed for second-order systems')
PSD=obj.filterPSD;

    Mz=PSD.Mz;
    Cz=PSD.Cz;
    Kz=PSD.Kz;

m=length(Mz);
    G= sparse(obj.n,m); G(obj.n,1)=1;


    M = [Mz, sparse(m,obj.n); sparse(obj.n,m), obj.M];
    C = [Cz, sparse(m,obj.n); sparse(obj.n,m), obj.C];
    K = [Kz, sparse(m,obj.n); -G, obj.K];

x=q(m+1:end);
xd=qd(m+1:end);

F_elastic = K * q + [sparse(m,1); obj.compute_fnl(x,xd)];

S=PSD.S;sigma=sqrt(S*2*pi); detT=obj.sde.timespan/obj.sde.num_points;
dF=randn*sqrt(detT)*sigma;
F_external = zeros(m+obj.n,1); F_external(m) =dF;
% F_external = obj.compute_fstochastic(t);

F_inertial = M * qdd;
F_damping = C * qd;

r = F_inertial + F_damping + F_elastic - F_external ;
drdqdd = M;
% initialization
dfdq=zeros(m+obj.n,m+obj.n); dfdqd=dfdq;
% compute the jacobian
dfdqd(m+1:end,m+1:end)=obj.compute_dfnldxd(x,xd);

dfdq(m+1:end,m+1:end)=obj.compute_dfnldx(x,xd);

drdqd = C + dfdqd;
drdq = K + dfdq;
%% 
% We use the following measure to comapre the norm of the residual $\mathbf{r}$
% 
% $$\texttt{c0} = \|\mathbf{M}\ddot{\mathbf{q}}\| + \|\mathbf{C}\dot{\mathbf{q}}\| 
% + \|\mathbf{F}(\mathbf{q})\| + \|\mathbf{F}_{ext}(t)\|$$
c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external);
end