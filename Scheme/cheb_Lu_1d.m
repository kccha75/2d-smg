% Performs multiplication L*u using Chebyshev collocation methods
% L is the spectral operator for the given PDE in the form
% u_xx+au_x+bu using
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
%
% Ouputs:
% A - L*u
%

function Lu=cheb_Lu_1d(v,pde,domain)

a=pde.a;
b=pde.b;

Lu=ifct(chebdiff(fct(v),2))+a.*ifct(chebdiff(fct(v),1))+b.*v;
Lu(1)=v(1);
Lu(end)=v(end);

end