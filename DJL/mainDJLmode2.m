clear;close all;%clc

% -------------------------------------------------------------------------
% DJL parameters
% -------------------------------------------------------------------------

epsilon=1;
alpha=epsilon^2;
mu=sqrt(epsilon);
L=1; % non-dimensionalised length scale of topography
u=0.1115;
KAI=20;
mode=2;

% N^2 function
N2=@(psi) psi;

% (N^2)'
N2d=@(psi) 0*psi+1;

% Linear part of N^2 (if exist)
lin=@(z,u) z/u^2;

% Nonlinear part of N^2
nonlin=@(z,v,u) -v.^2/u^2;

DJL.epsilon = epsilon;
DJL.alpha = alpha;
DJL.mu = mu;
DJL.L  = L;
DJL.u = u;
DJL.KAI = KAI;
DJL.mode=mode;

DJL.N2=N2;
DJL.N2d=N2d;
DJL.lin=lin;
DJL.nonlin=nonlin;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=DJLinitialise();

% Initial guess
v0=DJLv0(DJL,domain);

% Initialise PDE
[pde,domain]=DJLpdeinitialise(DJL,domain);

% Newton solve
[v,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);

% Continuation
[V,U]=naturalparametercontinuation(v,u,DJL,domain,cont_option);

dt=toc(time);
fprintf('Elapsed Time is %f s\n',dt)
% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------
X2=domain.X{1}/pi*KAI*L/mu;
Y2=(domain.X{2}+1)/2;

% Calculate momentum
P=trapI(V.^2,domain.dx{1}); % Integrate x
P=permute(P,[2,1,3]);
P=clenshaw_curtis(2*P/pi*KAI*L/mu); % Integrate y
P=permute(P,[3,1,2]);

plot(U,P)
xlabel('u');ylabel('Momentum');title('mode 1 DJL')

figure
contour(X2,Y2,Y2-V(:,:,end),100)
title("C=" + U(end))

figure;
plot(X2,Y2-V(:,:,end))
title("C=" + U(end))

% check dv<1 requirement
dv=2*ifct(chebdiff(fct(V(:,:,end)'),1));
max(dv(:))
min(dv(:))

fprintf('Domain is %d\n',KAI*L/mu)