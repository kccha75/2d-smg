clear;close all;%clc

% -------------------------------------------------------------------------
% fKdV parameters
% -------------------------------------------------------------------------

L=30;
d=3;
gamma=-0.5;
delta=1;
topography=@(x) sech(x).^2;

fKdV.L=L;
fKdV.d=d;
fKdV.gamma=gamma;
fKdV.delta=delta;
fKdV.topography=topography;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=fKdVinitialise();

ds=cont_option.ds;

% Initialise PDE
[pde,domain]=fKdVpdeinitialise(fKdV,domain);

v=0*domain.X; % initial guess

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
[pde,domain]=fKdVpdeinitialise(fKdV,domain); % update parameter

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
fprintf('Elapsed Time is %f s\n',dt)