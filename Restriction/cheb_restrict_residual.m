% Chebyshev Restriction operator in the x direction only (first component)
% assuming fine grid is 2x coarse grid points with 
%
% Adjoint operator adjustment applied
%
% Use for residual and RHS restriction 
% Boundary conditions applied (injection of fine grid)
%
% Inputs:
% vf - fine grid (Nx * Ny * ... size)
%
% Outputs:
% vc - coarse grid (ceil(Nx/2+1) * Ny * ... size)

function vc=cheb_restrict_residual(vf)

N=size(vf);

% Chebyshev transform
v_hat=fct(vf);

% For adjoint operator adjustment
v_hat(ceil(N(1)/2),:)=v_hat(ceil(N(1)/2),:)/2;

% Remove high frequency components and invert
vc=ifct(v_hat(1:ceil(N(1)/2),:));

% Apply BCs
vc(1,:)=vf(1,:);
vc(end,:)=vf(end,:);
% vc(1,:)=0;
% vc(end,:)=0;

% Reshape to correct size (required for 3d or higher)
N(1)=ceil(N(1)/2);
vc=reshape(vc,N);

end
