clear;%close all;%clc

% -------------------------------------------------------------------------
% DJL parameters PICK alpha / mu
% -------------------------------------------------------------------------
% fKdV solution type:
% 0 - 2sech^2 solution
% 1 - fKdV continuation plot!
DJL.soltype=3; 

mode=1; % mode solution
delta_star=1.5;%alpha=0.01; % topography height
gamma_star=0.1;% mu=0.7;
mu=0.90; % topography width scale
KAI=25; % fKdV domain, since L=200

% N^2 function
N2=@(psi) 1.0/(exp(1.0)-1)*exp(1.0*psi);N2=@(psi) 2*(-psi+1);

% (N^2)'
N2d=@(psi) 1.0/(exp(1.0)-1)*exp(1.0*psi);N2d=@(psi) -2+0*psi;

DJL.delta_star=delta_star;
DJL.gamma_star=gamma_star;
DJL.mu=mu;
DJL.mode=mode;
DJL.N2=N2;
DJL.N2d=N2d;
DJL.topography=@(X) sech(X+12).^2+sech(X-12).^2; % in KAI domain ...

DJL.KAI=KAI;

% -------------------------------------------------------------------------
time=tic;

% Initialise
[domain,option,cont_option]=DJLinitialise_topography();

% Initial guess
[DJL,fKdV,pdefkdv,domainfkdv,optionfkdv]=DJLv0_2topography(DJL,domain,option);
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

% Interpolate, leave NaN values for extrapolated
v0=interp2(H*(domain.X{2}+1)/2,domain.X{1},v0,YY,XX,'makima',NaN);

% Set BC here ... (but causes huge residual due to discont)
v0(:,end)=DJL.alpha*DJL.topography(domain.XX(:,end)*KAI/pi);

% Interpolate missing data (negative bump generally)
v0(isnan(v0))=griddata(YY(~isnan(v0)),XX(~isnan(v0)),v0(~isnan(v0)),YY(isnan(v0)),XX(isnan(v0)),'cubic');

% Fill missing (extrapolated) data (on boundary generally)
v0=fillmissing(v0,'makima',1);

% boundary layer for positive bump (smooths solution)
if DJL.alpha>0

    v0(2:end-1,end-round(0.05*size(domain.x{2})):end-1)=NaN;
    v0(isnan(v0))=griddata(YY(~isnan(v0)),XX(~isnan(v0)),v0(~isnan(v0)),YY(isnan(v0)),XX(isnan(v0)),'cubic');

end

% Initial residual
r=pde.f-(Lu(v0,pde,domain)+N2((domain.X{2}+1)/2-v0).*v0/DJL.u^2);
disp(rms(rms(r)))

% -------------------------------------------------------------------------
% Newton solve solution 1
% -------------------------------------------------------------------------

u1=DJL.u;
[v1,i,flag]=NewtonSolve(v0,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Newton solve solution 2 negative direction
% -------------------------------------------------------------------------

ds=1e-5;

DJL.u=u1-ds;u2=u1-ds;
[v2,i,flag]=NewtonSolve(v1,DJL,pde,domain,option);

if flag ==0

    fprintf('Initial Newton did not converge ...\n')
    return

end

% -------------------------------------------------------------------------
% Find Tabletop solution using Secant method
% -------------------------------------------------------------------------

[v,u,y,i,secantflag]=DJLtabletopsecant(v1,v2,u1,u2,DJL,pde,domain,option);

% % -------------------------------------------------------------------------
% % Adjust alpha for second solution
% % -------------------------------------------------------------------------
% 
% DJL.alpha=DJL.alpha+cont_option.ds;
% DJL.u=u;
% DJL.v=v;
% 
% % Conformal mapping and interpolation
% [DJL,domain]=conformalmapping(DJL,domain,option);
% 
% % Initialise PDE
% [DJL,pde,domain]=DJLpdeinitialise_topography(DJL,domain);
% 
% % Newton 1 here ...
% [v1,i,newtonflag1]=NewtonSolve(DJL.v,DJL,pde,domain,option);
%     
% % Newton 2 small delta adjustment
% DJL.u=DJL.u-1e-4;
% [v2,~,newtonflag2]=NewtonSolve(v1,DJL,pde,domain,option);
% 
% % -------------------------------------------------------------------------
% % Find Tabletop solution using Secant method
% % -------------------------------------------------------------------------
% 
% [vv,uu,y,secantflag2]=DJLtabletopsecant(v1,v2,u1,u2,DJL,pde,domain,option);
% 
% % Find dv
% dv=vv-v;
% dw=uu-u;
tic
% Continuation
if secantflag==1
    [V,U,W]=naturalparametercontalphaDJL(v,u,DJL.alpha,DJL,domain,option,cont_option);
end
toc
% % -------------------------------------------------------------------------
% % Newton solve solution 1
% % -------------------------------------------------------------------------
% 
% [v1,i,flag]=NewtonSolve(v,DJL,pde,domain,option);
% 
% if flag ==0
% 
%     fprintf('Initial Newton did not converge ...\n')
%     return
% 
% end