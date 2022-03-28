% Function calculates area approximation using trapezoidal rule
% in the x direction only (first component)
%
% Inputs:
% f - f(x) (Nx * Ny * ... size)
% x - vector of x
%
% Outputs:
% I (scalar in x direction)

function I=trapint(f,x)

N=size(f);

% Make column vector
x=x(:);

f=f(:,:);

% Formula for trapezoidal rule with summation through matrix multiplcation
I = 0.5*(x(2:end)-x(1:end-1))'*(f(2:end,:)+f(1:end-1,:));

% Scalar in x direction
N(1)=1;

% Reshape to correct size
I=reshape(I,N);

end