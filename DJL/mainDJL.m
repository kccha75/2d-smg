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

DJL.epsilon = epsilon;
DJL.alpha = alpha;
DJL.mu = mu;
DJL.L  = L;
DJL.u = u;
DJL.KAI = 20;
DJL.mode=1;
% -------------------------------------------------------------------------

[domain,option]=DJLinitialise();

v0=DJLv0(DJL,domain);

[pde,domain,option]=DJL_pde_initialise(DJL,domain,option);

[v,i,flag]=NewtonSolve(v0,pde,domain,option);

% continuation options?
cont_option.ds=0.001;
cont_option.ds_min=1e-6;
cont_option.ds_max=0.05;
cont_option.maxit= 8;
cont_option.N_opt=4;
cont_option.Newtontol=1e-10;
cont_option.steps=200;
cont_option.jacobian=@jacobian_DJL;

cont_option.operator=option.operator;
cont_option.preconditioner=option.preconditioner;
cont_option.prenumit=option.prenumit;
cont_option.tol=option.tol;

testDJLcontinuation
% DJLcontinuation

X2=domain.X{1}/pi*KAI*L/mu;
Y2=(domain.X{2}+1)/2;

plot(U,V)
xlabel('u');ylabel('Momentum');title('mode 1 DJL')

figure
contour(X2,Y2,Y2-v,100)
title("C=" + U(end))

figure;
plot(X2,Y2-v)
title("C=" + U(end))

% check dv<1 requirement
dv=2*ifct(chebdiff(fct(v'),1));
max(dv(:))
min(dv(:))

fprintf('Domain is %d\n',KAI*L/mu)

surf(v0)
min(v0(:))



