% Function finds fKdV solution for single sech^2 bump
%
% Inputs:
% DJL structure
% gamma - initial gamma value for continuation
% delta - delta in fKdV
%
% Ouputs:
% V - solution vector at each continuation steps
% U - parameter value at continuation steps
% fKdV - fKdV structure from DJL structure
% pde - fKdV pde parameters
% domain - fKdV domain
% option - fKdV solver options

function [V,U,fKdV,pde,domain,option]=fkdvsol(DJL,gamma,delta)

d=3; % nonlinear u^2 coefficient
gamma=min(gamma,-5); % initial gamma0

fKdV.L = 2*DJL.KAI;
fKdV.d = d;
fKdV.gamma = gamma;
fKdV.delta = delta;
fKdV.topography = DJL.topography;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=fKdVinitialise();

% Initialise PDE
[fKdV,pde,domain]=fKdVpdeinitialise(fKdV,domain);

% Initial guess
v=-0.75*sech(fKdV.L/(2*pi)*domain.X).^2; 

% -------------------------------------------------------------------------
% Newton solve for initial solution
% -------------------------------------------------------------------------
[v0,i,flag]=NewtonSolve(v,fKdV,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Newton solve for second solution
% -------------------------------------------------------------------------
ds=cont_option.ds;

fKdV.gamma=gamma+ds;

[fKdV,pde,domain]=fKdVpdeinitialise(fKdV,domain); % update parameter

[v1,i,flag]=NewtonSolve(v0,fKdV,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end
% -------------------------------------------------------------------------
% Set continuation initial conditions
v=v0;
u=gamma;
dv=(v1-v0)/ds;
du=1;

% Normalise
mag=sqrt(dot(dv,dv)+du^2);
dv=dv/mag;
du=du/mag;

% Continuation
[V,U]=pseudocontgamma(v,dv,u,du,fKdV,domain,option,cont_option);

dt=toc(time);
fprintf('fKdV Continuation Elapsed Time is %f s\n',dt)

end