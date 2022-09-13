clear;%close all;%clc

% -------------------------------------------------------------------------
% DJL parameters
% -------------------------------------------------------------------------
% fKdV solution type:
% 0 - 2sech^2 solution
% 1 - fKdV continuation plot!
DJL.soltype=1; 

% delta=0.01;
mode=1;
alpha=0.01;alpha=0;
mu=0.5;mu=1;
KAI=30;KAI=15;KAI=41;

% N^2 function
N2=@(psi) sech((psi-1)/1).^2;%N2=@(psi) psi;

% (N^2)'
N2d=@(psi) -2*sech((psi-1)/1).^2.*tanh((psi-1)/1);%N2d=@(psi) 1+0*psi;

DJL.alpha=alpha;
DJL.mode=mode;
DJL.N2=N2;
DJL.N2d=N2d;
% DJL.delta=delta;
DJL.topography=@(X) sech(X).^2; % in KAI domain ...

DJL.mu = mu;
DJL.KAI=KAI;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=DJLinitialise_topography();

% Initial guess
DJL=DJLv0_topography_test4(DJL,domain);

% Conformal mapping and interpolation
[DJL,domain]=conformalmapping(DJL,domain,option);

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

% set BC again here ... (but causes huge residual due to discont)
v0(:,end)=DJL.alpha*DJL.topography(domain.XX(:,end)*KAI/pi);

% initial residual ..?
r=pde.f-(Lu_2d(v0,pde,domain)+N2((domain.X{2}+1)/2-v0).*v0/DJL.u^2);
disp(rms(rms(r)))

% -------------------------------------------------------------------------
% Newton solve solution 1
% -------------------------------------------------------------------------
[v1,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end
% -------------------------------------------------------------------------
% Newton solve solution 2
% -------------------------------------------------------------------------
ds=cont_option.ds;
DJL.u=DJL.u+ds;
[v2,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end
% -------------------------------------------------------------------------

v=v1;
u=DJL.u-ds;
dv=(v2-v1)/ds;
du=1;

% Continuation
% [V,U]=naturalparametercontinuation(v,DJL.u,DJL,domain,cont_option);
[V,U]=pseudocontDJL(v,dv,u,du,DJL,domain,option,cont_option);

dt=toc(time);
fprintf('Elapsed Time is %f s\n',dt)

% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------

X2=domain.X{1}/pi*KAI/mu^2;
Y2=(domain.X{2}+1)/2;

% Calculate momentum
P=trapI(V.^2,domain.dx{1}); % Integrate x
P=permute(P,[2,1,3]);
P=clenshaw_curtis(2*P/pi*KAI/mu^2); % Integrate y
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