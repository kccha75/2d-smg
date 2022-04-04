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
mode=1;

DJL.mode=1;
DJL.epsilon = epsilon;
DJL.alpha = alpha;
DJL.mu = mu;
DJL.L  = L;
DJL.u = u;

DJL.KAI = 20;
% -------------------------------------------------------------------------

[domain,option]=DJLinitialise();

v0=DJLv0(DJL,domain);

[pde,domain]=DJL_pde_initialise(DJL,domain);

v=DJLsolve(v0,pde,domain,option);
% 
% X2=X/pi*KAI*L/mu;
% Y2=(Y+1)/2;
% 
% DJLcontinuation


% surf(X2,Y2,v0);title('initial guess')
% figure;surf(X2,Y2,v);title('DJL solution');xlabel('x');ylabel('z')
% figure;surf(X2,Y2,v0-v);title('initial guess error');xlabel('x');ylabel('z')

fprintf('Domain is %d\n',KAI*L/mu)

surf(v0)
min(v0(:))