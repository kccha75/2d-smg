% Function calculates area approximation using trapezoidal rule
% in the x direction only (first component)
%
% Inputs:
% f - f(x) (Nx * Ny * ... size)
% dx - step size (assumed constant)
%
% Outputs:
% I (scalar in x direction)

function I=trapI(f,dx)

N=size(f);

f=f(:,:);

% Formula for trapezoidal rule with summation in x direction
I = sum(dx*f);

% Scalar in x direction
N(1)=1;

% Reshape to correct size
I=reshape(I,N);

end