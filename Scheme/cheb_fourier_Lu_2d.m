% Performs multiplication L*u using Cheb-Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% au_xx+bu_yy+cu=f
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% domain.k - wave number
%
% Ouputs:
% A - L*u
%
function Lu=cheb_fourier_Lu_2d(v,pde,domain)

kx=domain.k{1};
ky=domain.k{2};

a=pde.a;
b=pde.b;
c=pde.c;

% au_xx term
Lu1=ifct(chebdiff(fct(v),2));
% Lu1(1,:)=0;Lu1(end,:)=0;
% bu_yy term
Lu2=real(ifft(-ky.^2.*fft(v')));


Lu=a.*Lu1+b.*Lu2'+c.*v;

% Apply BCs
Lu(1,:)=v(1,:);
Lu(end,:)=v(end,:);

end