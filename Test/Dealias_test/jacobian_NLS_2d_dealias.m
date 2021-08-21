% Gives the Jacobian for the PDE in the form
% -(au_x)_x-b(u_y)_y+cu+du^3
%
% Jacobian is -(ae_x)_x-b(e_y)_y+(c+3dv^2)e
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

function jacobian=jacobian_NLS_2d_dealias(v,pde,domain)

jacobian.a=pde.a;
jacobian.b=pde.b;
jacobian.c=pde.c+3*pde.d.*dealias_2d(v,v);
jacobian.f=pde.f-(fourier_Lu_2d_dealias(v,pde,domain)+pde.d.*dealias_2d(dealias_2d(v,v),v));

end