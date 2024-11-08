% Continuation uses mu09continuationstarting results, pick a mu and do a
% normal continuation in alpha (does not use initial guess)

% load data here! 
clear; 
load('gammastar05/exp/U.mat') % initial alpha
alpha=U(1);

% change if want to go positive or negative direction
% load('mu09continuationstarting/U.mat') % mu
% load('mu09continuationstarting/V.mat') % solution
% load('mu09continuationstarting/W.mat') % alpha

load('gammastar05/expcontpos/U.mat')
load('gammastar05/expcontpos/V.mat')
load('gammastar05/expcontpos/W.mat')

% -------------------------------------------------------------------------
% DJL parameters PICK alpha / mu
% -------------------------------------------------------------------------
% fKdV solution type:
% 0 - 2sech^2 solution
% 1 - fKdV continuation plot!
DJL.soltype=3; 

mu=1.0; % topography width scale
KAI=25; % fKdV domain, since L=200

% N^2 function
N2=@(psi) exp(-1.0*psi)/((exp(-1.0)-1)/-1.0);

% (N^2)'
N2d=@(psi) -1.0*exp(-1.0*psi)/((exp(-1.0)-1)/-1.0);

% find mu index or approximately ...
[~,muindex]=min(abs(U-mu));

DJL.mu=mu;
DJL.N2=N2;
DJL.N2d=N2d;
DJL.topography=@(X) sech(X+12).^2+sech(X-12).^2; % in KAI domain ...
DJL.alpha=alpha; % continuation point!
DJL.KAI=KAI;
DJL.Lx=2*KAI/mu^2;
% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=DJLinitialise_topography();

% Conformal mapping and interpolation
[DJL,domain]=conformalmapping(DJL,domain,option);

% Length scales in DJL coordinates
Lx=DJL.Lx;
Ly=DJL.Ly;

XX=domain.XX;
YY=domain.YY;
jac=domain.jac;
H=domain.H;

% Initialise PDE
[DJL,pde,domain]=DJLpdeinitialise_topography(DJL,domain);

v=V(:,:,muindex); % continuation point!
u=W(muindex); % continuation point!
DJL.u=u;
% -------------------------------------------------------------------------
% Newton solve solution 1
% -------------------------------------------------------------------------

u1=DJL.u;
[v1,i,flag]=NewtonSolve(v,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Newton solve solution 2 negative direction
% -------------------------------------------------------------------------

ds=1e-5;

DJL.u=u1-ds;u2=u1-ds;
[v2,i,flag]=NewtonSolve(v1,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Find Tabletop solution using Secant method
% -------------------------------------------------------------------------

[v,u,y,i,secantflag]=DJLtabletopsecant(v1,v2,u1,u2,DJL,pde,domain,option);

% continuation
if secantflag==1
    [V,U,W]=naturalparametercontalphaDJL(v,u,DJL.alpha,DJL,domain,option,cont_option);
end
