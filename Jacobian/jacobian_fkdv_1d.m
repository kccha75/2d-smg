% Gives the Jacobian for the PDE in the form
% (-u_xx+au_x+bu+cu^2)
%
% Jacobian is (-e_xx+ae_x+(b+2cv)e)
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

function jacobian=jacobian_fkdv_1d(v,pde,domain)

jacobian.a=pde.a;
jacobian.b=pde.b+2*pde.c.*v;
jacobian.c=0;
jacobian.f=pde.f-(fourier_Lu_1d(v,pde,domain)+pde.c.*v.^2);

end