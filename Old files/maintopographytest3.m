clear;%close all;%clc

% -------------------------------------------------------------------------
% DJL parameters PICK delta / mu
% -------------------------------------------------------------------------
% fKdV solution type:
% 0 - 2sech^2 solution
% 1 - fKdV continuation plot!
DJL.soltype=1; 
delta=0.01; % DJL perturbation
% try delta_star=3
% gamma_star=-0.5, solution 3 :)
mode=1; % mode solution
mu=0.7; % topography width scale
KAI=30;KAI=15; % fKdV domain

% N^2 function
N2=@(psi) sech((psi-1)/1).^2;%N2=@(psi) psi;

% (N^2)'
N2d=@(psi) -2*sech((psi-1)/1).^2.*tanh((psi-1)/1);%N2d=@(psi) 1+0*psi;

DJL.mode=mode;
DJL.N2=N2;
DJL.N2d=N2d;
DJL.delta=delta;
DJL.topography=@(X) sech(X).^2; % in KAI domain ...

DJL.mu = mu;
DJL.KAI=KAI;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=DJLinitialise_topography();

% Initial guess
DJL=DJLv0_topography_test3(DJL,domain);

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

% Interpolate
v0=interp2(H*(domain.X{2}+1)/2,domain.X{1},DJL.v,YY,XX,'spline');

% Set BC again here ... (but causes huge residual due to discont)
v0(:,end)=DJL.alpha*DJL.topography(domain.XX(:,end)*KAI/pi);

% Initial residual ..?
r=pde.f-(Lu_2d(v0,pde,domain)+N2((domain.X{2}+1)/2-v0).*v0/DJL.u^2);
disp(rms(rms(r)))

% Newton solve
[v,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% Continuation
[V,U]=naturalparametercontinuation(v,DJL.u,DJL,domain,option,cont_option);

dt=toc(time);
fprintf('Elapsed Time is %f s\n',dt)

% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------

X2=domain.X{1}/pi*KAI;
Y2=(domain.X{2}+1)/2;

% Calculate momentum
P=trapI(V.^2,domain.dx{1}); % Integrate x
P=permute(P,[2,1,3]);
P=clenshaw_curtis(2*P/pi*KAI); % Integrate y
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

fprintf('Domain is %d\n',KAI/mu^2)