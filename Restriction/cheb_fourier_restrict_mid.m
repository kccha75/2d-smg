% 2D Cheb-Fourier Restriction operator assuming fine grid is 2x
% coarse grid points
%
% Assumes using midpoint method with adjustment
%
% Inputs:
% vf - fine grid (Nx * Ny size)
%
% Outputs:
% vc - fine grid restriction (Nx/2 * Ny/2 size)

function vc=cheb_fourier_restrict_mid(vf)

% Chebyshev restrict in x
vc=cheb_restrict(vf);

% Fourier restrict in y
vc=fourier_restrict_mid(vc');
vc=vc';

end