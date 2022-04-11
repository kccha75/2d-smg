% Gives the Jacobian for the DJL PDE in the form
%
% A_xx+A_zz+N^2(z-A)*A/u^2=0
%
% Jacobian is e_xx+e_zz+(eN^2(z-v)-e*v*dN^2(z-v))/u^2
%
% Inputs:
% DJL.N2 - N^2 function
% DJL.N2d - dN^2 function
% pde.a - coefficient in front of A_xx
% pde.b - coefficient in front of A_zz
% pde.c - coefficient in front of A
% pde.u - wave speed
% domain.z
%
% Ouputs:
% jacobian.a
% jacobian.b
% jacobian.c - coefficient (N^2(z-v)-v*dN^2(z-v))/u^2
% jacobian.f - nonlinear residual

function J=jacobianDJL(v,DJL,pde,domain)

N2=DJL.N2;
N2d=DJL.N2d;

u=pde.u;

z=(domain.X{2}+1)/2;

J.a=pde.a;
J.b=pde.b;
J.c=(N2(z-v)-v.*N2d(z-v))/u^2;
J.f=pde.f-(Lu_2d(v,pde,domain)+N2(z-v).*v/u^2); % nonlinear residual

end