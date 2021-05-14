% 2D Fourier Restriction operator assuming fine grid is 2x
% coarse grid points
%
% Assumes using midpoint method with adjustment
%
% Inputs:
% vf - fine grid (Nx * Ny size)
%
% Outputs:
% vc - fine grid restriction (Nx/2 * Ny/2 size)

function vc=fourier_restrict_2d_mid(vf)

% FFT + shift in x
vc=fourier_restrict_1d_mid(vf);

% FFT + shift in y
vc=fourier_restrict_1d_mid(vc');

vc=vc';

end
