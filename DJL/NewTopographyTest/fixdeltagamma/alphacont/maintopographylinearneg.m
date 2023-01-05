clear;%close all;%clc

% -------------------------------------------------------------------------
% DJL parameters PICK alpha / mu
% -------------------------------------------------------------------------
% fKdV solution type:
% 0 - 2sech^2 solution
% 1 - fKdV continuation plot!
DJL.soltype=1; 

mode=1; % mode solution
delta_star=1.5;%alpha=0.01; % topography height
gamma_star=0.25;% mu=0.7;
mu=0.80; % topography width scale
KAI=30;KAI=20; % fKdV domain

% N^2 function
N2=@(psi) 2*(-psi+1);%N2=@(psi) psi;

% (N^2)'
N2d=@(psi) -2+0*psi;%N2d=@(psi) 1+0*psi;

DJL.delta_star=delta_star;
DJL.gamma_star=gamma_star;
DJL.mu=mu;
DJL.mode=mode;
DJL.N2=N2;
DJL.N2d=N2d;
DJL.topography=@(X) sech(X).^2; % in KAI domain ...

DJL.KAI=KAI;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=DJLinitialise_topography();

% Initial guess
[DJL,fKdV,pdefkdv,domainfkdv,optionfkdv]=DJLv0_topography(DJL,domain,option);
v0=DJL.v;

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
v0=interp2(H*(domain.X{2}+1)/2,domain.X{1},v0,YY,XX,'spline');

% Set BC again here ... (but causes huge residual due to discont)
v0(:,end)=DJL.alpha*DJL.topography(domain.XX(:,end)*KAI/pi);

% Initial residual
r=pde.f-(Lu(v0,pde,domain)+N2((domain.X{2}+1)/2-v0).*v0/DJL.u^2);
disp(rms(rms(r)))

% -------------------------------------------------------------------------
% Newton solve solution 1
% -------------------------------------------------------------------------

u1=DJL.alpha;
[v1,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Newton solve solution 2 negative direction
% -------------------------------------------------------------------------

ds=cont_option.ds;ds=0.01;

DJL.alpha=u1-ds;

% Mapping
[DJL,domain]=conformalmapping(DJL,domain,option);
[DJL,pde,domain]=DJLpdeinitialise_topography(DJL,domain);

[v2,i,flag]=NewtonSolve(v1,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Continuation DJL negative direction
% -------------------------------------------------------------------------

v=v1;
u=u1;
dv=(v2-v1)/ds;
du=-1; % +1 for positive direction, -1 for negative direction

% Continuation
[Vneg,Uneg]=pseudocontDJL2(v,dv,u,du,DJL,domain,option,cont_option);

dt=toc(time);
fprintf('Elapsed Time is %f s\n',dt)

% -------------------------------------------------------------------------
% Newton solve solution 2 positive direction
% -------------------------------------------------------------------------

DJL.u=u1+ds;
[v2,i,flag]=NewtonSolve(v1,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Continuation DJL positive direction
% -------------------------------------------------------------------------

v=v1;
u=u1;
dv=(v2-v1)/ds;
du=1; % +1 for positive direction, -1 for negative direction

% Continuation
[Vpos,Upos]=pseudocontDJL2(v,dv,u,du,DJL,domain,option,cont_option);

dt=toc(time);
fprintf('Elapsed Time is %f s\n',dt)


% -------------------------------------------------------------------------
% PLOTS
% -------------------------------------------------------------------------
U=[fliplr(Uneg) Upos];
V=cat(3,flip(Vneg,3),Vpos);

X2=domain.X{1}/pi*KAI/mu;
Y2=(domain.X{2}+1)/2;

% Calculate momentum
P=trapI(V.^2,2*KAI/mu/domain.N(1)); % Integrate x
P=permute(P,[2,1,3]);
P=clenshaw_curtis(P)/2; % Integrate y
P=permute(P,[3,1,2]);

% u vs momentum
plot(U,P)
xlabel('u');ylabel('Momentum');title('mode 1 DJL')

% Contour of final solution
figure
contour(X2,Y2,Y2-Vneg(:,:,end),100)
title("Negative C=" + Uneg(end))

figure
contour(X2,Y2,Y2-Vpos(:,:,end),100)
title("Positive C=" + Upos(end))

% Plot(s) of final solution
figure;
plot(X2,Y2-Vneg(:,:,end))
title("Negative C=" + Uneg(end))

% Plot(s) of final solution
figure;
plot(X2,Y2-Vpos(:,:,end))
title("Positive C=" + Upos(end))

% check dv<1 requirement
dv=2*real(ifct(chebdiff(real(fct(transpose(Vneg(:,:,end)))),1)));
max(dv(:))
min(dv(:))
dv=2*real(ifct(chebdiff(real(fct(transpose(Vpos(:,:,end)))),1)));
max(dv(:))
min(dv(:))

fprintf('Domain is %d\n',KAI/mu^2)

% -------------------------------------------------------------------------
% KdV solution 1
% -------------------------------------------------------------------------

B1=fKdV.B;
delta=fKdV.delta;

% -------------------------------------------------------------------------
% KdV solve 2 negative direction
% -------------------------------------------------------------------------

fKdV.delta=delta-ds;

[fKdV,pdefkdv,domainfkdv]=fKdVpdeinitialise(fKdV,domainfkdv); % update parameter

[B2,i,flag]=NewtonSolve(B1,fKdV,pdefkdv,domainfkdv,optionfkdv);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Continuation fKdV negative direction
% -------------------------------------------------------------------------

B=B1;
D=delta;
dB=(B2-B1)/ds;
dD=-1; % +1 for positive direction, -1 for negative direction

[Bneg,Dneg]=pseudocontdelta(B,dB,D,dD,fKdV,domainfkdv,optionfkdv,cont_option);

% -------------------------------------------------------------------------
% Newton solve solution 2 positive direction
% -------------------------------------------------------------------------

fKdV.delta=delta+ds;

[fKdV,pdefkdv,domainfkdv]=fKdVpdeinitialise(fKdV,domainfkdv); % update parameter

[B2,i,flag]=NewtonSolve(B1,fKdV,pdefkdv,domainfkdv,optionfkdv);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Continuation fKdV positive direction
% -------------------------------------------------------------------------

B=B1;
D=delta;
dB=(B2-B1)/ds;
dD=1; % +1 for positive direction, -1 for negative direction

[Bpos,Dpos]=pseudocontdelta(B,dB,D,dD,fKdV,domainfkdv,optionfkdv,cont_option);

% -------------------------------------------------------------------------
% Combine!
% -------------------------------------------------------------------------

D=[fliplr(Dneg) Dpos];
B=cat(2,flip(Bneg,2),Bpos);


% -------------------------------------------------------------------------
% 0th order DJL approx
% -------------------------------------------------------------------------

% Find A
A=6*DJL.s*mu^2/DJL.r*B;

A=reshape(A,[size(A,1),1,size(A,2)]);

% Find zeta from fkdv
zeta0=pagemtimes(A,DJL.phi');

% Conformal map!
for jj=1:size(zeta0,3)
    zeta0(:,:,jj)=interp2(H*(domain.X{2}+1)/2,domain.X{1},zeta0(:,:,jj),YY,XX,'spline');
    zeta0(:,end,jj)=DJL.alpha*DJL.topography(domain.XX(:,end)*KAI/pi);
end
% -------------------------------------------------------------------------
% 1st order DJL approx
% -------------------------------------------------------------------------

% A_xx in x domain (KAI/mu)
A_xx=ifft(-(pi/(1/mu*KAI)*domain.k{1}).^2.*fft(A));

% an
an=pagemtimes(A_xx,DJL.a1)+pagemtimes(A.^2,DJL.a2)+pagemtimes(DJL.alpha*DJL.b,DJL.a3);

% n=N case
an(:,mode,:)=0;

zeta1=pagemtimes(an,DJL.phis');

% Back to zai coordinates
zai=zeta1+DJL.alpha*DJL.b*(1-(domain.x{2}+1)/2)';

zeta=zeta0+zai;

% Conformal map!
for jj=1:size(zeta,3)
    zeta(:,:,jj)=interp2(H*(domain.X{2}+1)/2,domain.X{1},zeta(:,:,jj),YY,XX,'spline');
    zeta(:,end,jj)=DJL.alpha*DJL.topography(domain.XX(:,end)*KAI/pi);
end

% -------------------------------------------------------------------------
% Plots!
% -------------------------------------------------------------------------

% Calculate momentum 0th order
P2=trapI(zeta0.^2,2*KAI/mu/domain.N(1)); % Integrate x
P2=permute(P2,[2,1,3]);
P2=clenshaw_curtis(P2)/2; % Integrate y
P2=permute(P2,[3,1,2]);

% Calculate momentum 1st order
P3=trapI(zeta.^2,2*KAI/mu/domain.N(1)); % Integrate x
P3=permute(P3,[2,1,3]);
P3=clenshaw_curtis(P3)/2; % Integrate y
P3=permute(P3,[3,1,2]);

% u vs momentum
figure;
plot(DJL.C+D*DJL.s*mu^2,P2)
xlabel('delta');ylabel('Momentum');title('fKdV mode 1')

% % u vs momentum
% figure;
% plot(DJL.C+D*DJL.s*mu^2,P3)
% xlabel('delta');ylabel('Momentum');title('fKdV mode 1')

% Compare all 
figure;
plot(DJL.C+D*DJL.s*mu^2,P2,DJL.C+D*DJL.s*mu^2,P3,U,P)
xlabel('delta');ylabel('Momentum');title('DJL momentums')
legend('0th order fKdV approx of DJL','1st order fKdV approx of DJL','Exact DJL')

% figure;
% plot(DJL.C+D*DJL.s*mu^2,P2,U,P)
% xlabel('delta');ylabel('Momentum');title('DJL momentums')
% legend('0th order fKdV approx of DJL','Exact DJL')