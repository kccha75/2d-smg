% Chebyshev Prolongation operator on x direction (first component)
% assuming fine grid is 2x coarse grid points in x
%
% Inputs:
% vc - coarse grid (Nx * Ny * ... size)
%
% Outputs:
% vc - coarse grid prolongation (Nx*2 * Ny * ... size)

function vf=cheb_prolong(vc)

N=size(vc);

% Chebyshev transform
vc_hat=fct(vc);

% Filter high frequency and invert
vf=ifct([vc_hat;zeros(N)]);

end
