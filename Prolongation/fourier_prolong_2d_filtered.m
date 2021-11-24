% 2D Fourier Prolongation operator with x and y sweeps
% assuming fine grid is 2x coarse grid points
%
% Assumed filtered for highest mode
%
% Inputs:
% vc - coarse grid (Nx * Ny size)
%
% Outputs:
% vc - coarse grid prolongation (Nx*2 * Ny*2 size)

function vf=fourier_prolong_2d_filtered(vc)

% FFT + shift in x
vf=fourier_prolong_filtered(vc);

% FFT + shift in y
vf=fourier_prolong_filtered(vf');

vf=vf';

end