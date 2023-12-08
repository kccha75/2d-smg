clear;%close all;%clc

% -------------------------------------------------------------------------
% DJL parameters
% -------------------------------------------------------------------------

epsilon=1;%epsilon=10;
alpha=epsilon^2; % no dependence on alpha (since no topography ...)
mu=sqrt(epsilon);% no dependence on mu at all for v0... but does have dependence on newton ...
u=0.33;
mode=1;

% N^2 function
N2=@(psi) 1/(exp(1)-1)*exp(psi);

% (N^2)'
N2d=@(psi) 1/(exp(1)-1)*exp(psi);

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
[v0,DJL]=DJLv0(DJL,domain);

% Initialise PDE
[pde,domain]=DJLpdeinitialise(DJL,domain);

disp(rms(rms(pde.f-(Lu(v0,pde,domain)+N2((domain.X{2}+1)/2-v0).*v0/DJL.u^2))))

% Newton solve
tic;[v,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);toc;

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% Continuation
[Vpos,Upos]=naturalparametercontinuation(v,u,DJL,domain,option,cont_option);

% Negative direction
cont_option.direction=-1;
[Vneg,Uneg]=naturalparametercontinuation(v,u,DJL,domain,option,cont_option);

% Combine
U=[fliplr(Uneg) Upos];
V=cat(3,flip(Vneg,3),Vpos);

dt=toc(time);
fprintf('Elapsed Time is %f s\n',dt)
% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------
KAI=DJL.KAI;

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