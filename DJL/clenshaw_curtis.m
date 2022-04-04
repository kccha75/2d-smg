% Clenshaw Curtis intergration in the x direction only (first component)
%
% Inputs:
% f - f(x) at Chebyshev points (Nx * Ny * ... size)
%
% Ouputs:
% I - Integral (scalar in x direction)

function I = clenshaw_curtis(f)

N=size(f);
M=N(1)-1;

f=f(:,:);

% FFT
g = real(fft(f([1:M+1 M:-1:2],:)/(2*M)));

% Chebyshev coefficients
a = [g(1,:); g(2:M,:)+g(2*M:-1:M+2,:); g(M+1,:)];

% Weight vector
w = zeros(M+1,1);
w(1:2:end) = 2./(1-(0:2:M).^2);

% Integral
I = w'*a;

% Scalar in x direction
N(1)=1;

% Reshape to correct size
I=reshape(I,N);