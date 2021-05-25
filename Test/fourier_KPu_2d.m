% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% -au_xx+cu+b*intint(u_yy)dxdx using
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% pde.d
% domain.k - wave number
% domain.N
%
% Ouputs:
% A - L*u
%

function Lu=fourier_KPu_2d(v,pde,domain)

kx=domain.k(:,1);
ky=domain.k(:,2);
kx(1)=1e-6;
[KX,KY]=ndgrid(kx,ky);

a=pde.a;
b=pde.b;
c=pde.c;

Kux=-a.*ifft(-kx.^2.*fft(v))+c.*v;

Kuy=b.*ifft2(KY.^2./KX.^2.*fft2(v));

Lu=Kux+Kuy;

end