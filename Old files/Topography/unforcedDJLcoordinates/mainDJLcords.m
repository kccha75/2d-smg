clear;close all;%clc

% -------------------------------------------------------------------------
% DJL parameters
% -------------------------------------------------------------------------

epsilon=0.5; % appears no epsilon dependence on gamma_hat_star, delta delta_hat_star ...
alpha=epsilon^2; % no dependence on alpha (since no topography ...)
mu=sqrt(epsilon); % appears no mu dependence on gamma_hat_star, delta delta_hat_star ...
u=0.24;
mode=1;

% N^2 function
N2=@(psi) psi;

% (N^2)'
N2d=@(psi) 0*psi+1;

DJL.epsilon = epsilon;
DJL.alpha = alpha;
DJL.mu = mu;
DJL.u = u;
DJL.mode=mode;

DJL.N2=N2;
DJL.N2d=N2d;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=DJLinitialise();

% Initial guess
[v0,DJL]=DJLv0cords(DJL,domain);

% Initialise PDE
[pde,domain]=DJLpdeinitialisecords(DJL,domain);

disp(rms(rms(pde.f-(Lu_2d(v0,pde,domain)+N2((domain.X{2}+1)/2-v0).*v0/DJL.u^2))))

% Newton solve
[v,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% Continuation
[V,U]=naturalparametercontinuation(v,u,DJL,domain,cont_option);

dt=toc(time);
fprintf('Elapsed Time is %f s\n',dt)
% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------
KAI=DJL.KAI;

X2=domain.X{1}/pi*KAI/mu;
Y2=(domain.X{2}+1)/2;

% Calculate momentum
P=trapI(V.^2,domain.dx{1}); % Integrate x
P=permute(P,[2,1,3]);
P=clenshaw_curtis(2*P/pi*KAI/mu); % Integrate y
P=permute(P,[3,1,2]);

% u vs momentum
plot(U,P)
xlabel('u');ylabel('Momentum');title('mode 1 DJL')

% Contour of final solution
figure
contour(X2,Y2,Y2-V(:,:,end),100)
title("C=" + U(end))

% Plot(s) of final solution
figure;
plot(X2,Y2-V(:,:,end))
title("C=" + U(end))

% check dv<1 requirement
dv=2*ifct(chebdiff(fct(V(:,:,end)'),1));
max(dv(:))
min(dv(:))

fprintf('Domain is %d\n',KAI/mu)