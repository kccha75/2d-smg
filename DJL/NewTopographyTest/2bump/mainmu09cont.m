% Continuation on a specific point on mu09 continuation. attempts to
% continuation negatively (from 0.9 towards -infinity) to overturning
% can pick positive direction, see ds

% load data here!
clear; 
% load('mu09results/U.mat')
% load('mu09results/V.mat')
% load('mu09results/W.mat')
% 
load('gammastar05/exp/U.mat')
load('gammastar05/exp/V.mat')
load('gammastar05/exp/W.mat')

% load('U.mat')
% load('V.mat')
% load('W.mat')

% continuation point in data
contpoint=1;

% -------------------------------------------------------------------------
% DJL parameters PICK alpha / mu
% -------------------------------------------------------------------------
% fKdV solution type:
% 0 - 2sech^2 solution
% 1 - fKdV continuation plot!
DJL.soltype=3; 

mu=0.7; % topography width scale
KAI=25; % fKdV domain, since L=200

% N^2 function
N2=@(psi) exp(-1.0*psi)/((exp(-1.0)-1)/-1.0);

% (N^2)'
N2d=@(psi) -1.0*exp(-1.0*psi)/((exp(-1.0)-1)/-1.0);

DJL.mu=mu;
DJL.N2=N2;
DJL.N2d=N2d;
DJL.topography=@(X) sech(X+12).^2+sech(X-12).^2; % in KAI domain ...
DJL.alpha=U(contpoint); % continuation point!
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

v=V(:,:,contpoint); % continuation point!
u=W(contpoint); % continuation point!

% pick +ve or -ve direction
cont_option.ds=-cont_option.ds;

[V,U,W]=naturalparametercontalphaDJLmu(v,u,mu,DJL,domain,option,cont_option);
