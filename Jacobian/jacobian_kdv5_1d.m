% Gives the Jacobian for the PDE in the form
% (-u_xxxx+au_xx+bu+cu^2)
%
% Jacobian is (-e_xxxx+ae_xx+2bve)
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

function jacobian=jacobian_kdv5_1d(v,pde,domain)

jacobian.a=pde.a;
jacobian.b=pde.b+2*pde.c.*v;
jacobian.c=0;
jacobian.f=pde.f-fourier_Lu_1d_kdv5(v,pde,domain);

end