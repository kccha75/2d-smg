function vc=fourier_cheb_restrict(vf)

[Nx,Ny]=size(vf);

% Fourier restrict in x
vc=fourier_restrict_1d(vf);

% Chebyshev restrict in y
vc=cheb_restrict_1d(vc');
vc=vc';

end