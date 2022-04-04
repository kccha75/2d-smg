% Gives the Jacobian for the DJL PDE in the form
%
% A_xx+A_zz+(z-A)*A/u^2=0
%
% Jacobian is e_xx+e_zz+e(z-2v)/u^2
%
% Inputs:
% pde.a - coefficient in front of A_xx
% pde.b - coefficient in front of A_zz
% pde.c - coefficient in front of A
% pde.u - wave speed
% domain
%
% Ouputs:
% jacobian.a
% jacobian.b
% jacobian.c
% jacobian.f - nonlinear residual

function jacobian=jacobian_DJL(v,pde,domain)

jacobian.a=pde.a;
jacobian.b=pde.b;
jacobian.c=pde.c-2*v/pde.u^2;
jacobian.f=pde.f-(Lu_2d(v,pde,domain)-v.^2/pde.u^2);

end