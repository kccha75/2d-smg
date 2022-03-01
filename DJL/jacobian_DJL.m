% Gives the Jacobian for the PDE in the form
%
% A_xx+A_zz+(z-A)*A/u^2=0
%
% Jacobian is e_xx+e_zz+(ze-2ev)/u^2
%
% Inputs:
% pde.a
% pde.b
% pde.c
% domain
%
% Ouputs:
% jacobian.a
% jacobian.b
% jacobian.c

function jacobian=jacobian_DJL(v,pde,domain)

jacobian.a=pde.a;
jacobian.b=pde.b;
jacobian.c=???;
jacobian.f=???;

end