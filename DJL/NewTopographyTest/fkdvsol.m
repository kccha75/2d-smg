function [V,U,fKdV,pde,domain,option]=fkdvsol(DJL,gamma,delta)

d=3; % nonlinear u^2 coefficient
gamma=min(gamma,-5); % initial gamma0

fKdV.L=2*DJL.KAI;
fKdV.d=d;
fKdV.gamma=gamma;
fKdV.delta=delta;
fKdV.topography=DJL.topography;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=fKdVinitialise();

ds=cont_option.ds;

% Initialise PDE
[fKdV,pde,domain]=fKdVpdeinitialise(fKdV,domain);

v=-0.75*sech(fKdV.L/(2*pi)*domain.X).^2; % initial guess

% -------------------------------------------------------------------------
% Newton solve for initial solution
[v0,i,flag]=NewtonSolve(v,fKdV,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Newton solve for second solution

fKdV.gamma=gamma+ds;
[fKdV,pde,domain]=fKdVpdeinitialise(fKdV,domain); % update parameter

[v1,i,flag]=NewtonSolve(v0,fKdV,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end
% -------------------------------------------------------------------------
v=v0;
u=gamma;
dv=(v1-v0)/ds;
du=1;

% Normalise
mag=sqrt(dot(dv,dv)+du^2);
dv=dv/mag;
du=du/mag;

% Continuation
[V,U]=pseudocont(v,dv,u,du,fKdV,domain,option,cont_option);

dt=toc(time);
fprintf('fKdV Continuation Elapsed Time is %f s\n',dt)

end