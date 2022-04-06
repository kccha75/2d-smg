clear;close all;%clc

% -------------------------------------------------------------------------
% DJL parameters
% -------------------------------------------------------------------------

epsilon=1;
alpha=epsilon^2;
mu=sqrt(epsilon);
L=1; % non-dimensionalised length scale of topography
u=0.24;
KAI=20;
mode=2;

DJL.epsilon = epsilon;
DJL.alpha = alpha;
DJL.mu = mu;
DJL.L  = L;
DJL.u = u;
DJL.KAI = KAI;
DJL.mode=mode;

% -------------------------------------------------------------------------
% Initialise
[domain,option,cont_option]=DJLinitialise();

% Initial guess
v0=DJLv0(DJL,domain);

% Initialise PDE
[pde,domain,option]=DJL_pde_initialise(DJL,domain,option);

% Newton solve
[v,i,flag]=NewtonSolve(v0,pde,domain,option);

% Continuation
[V,U]=naturalparametercontinuation(v,u,DJL,domain,cont_option);

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