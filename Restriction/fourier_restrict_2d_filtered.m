% 2D Fourier Restriction operator assuming fine grid is 2x
% coarse grid points
%
% Assumed filtered for highest mode
%
% Inputs:
% vf - fine grid (Nx * Ny size)
%
% Outputs:
% vc - fine grid restriction (Nx/2 * Ny/2 size)

function vc=fourier_restrict_2d_filtered(vf)

% FFT + shift in x
vc=fourier_restrict_filtered(vf);

% FFT + shift in y
vc=fourier_restrict_filtered(vc');

vc=vc';

end