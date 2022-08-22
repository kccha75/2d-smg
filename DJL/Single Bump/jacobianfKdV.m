% Gives the Jacobian for the ODE (quadratic nonlinearity) in the form
%
% au_xx+bu_x+cu+du^2=f
%
% Jacobian is ae_xx+be_x+(c+2dv)e=f-(av_xx+bv_x+cv+dv^2)
%
% Inputs:
% v - best estimate
% fKdV.d - coefficient d in front of u^2
% pde.a - coefficient in front of u_xx
% pde.b - coefficient in front of u_x
% pde.c - coefficient in front of u
%
% Ouputs:
% jacobian.a
% jacobian.b
% jacobian.c - coefficient c+2*d*v
% jacobian.f - nonlinear residual

function J=jacobianfKdV(v,fKdV,pde,domain)

d=fKdV.d;

J.a=pde.a;
J.b=pde.b;
J.c=pde.c+2*d.*v;
J.f=pde.f-(Lu(v,pde,domain)+d.*v.^2); % nonlinear residual

end