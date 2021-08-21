% Fast Chebyshev differentiation in the x direction only (first component)
%
% Note: Must be at Chebyshev collocation points
%
% Inputs:
% v - function (Nx * Ny * ... size)
%
% Ouputs:
% dv - derivative function (Nx * Ny * ... size)
%

function w = fchd(v)

v = v(:)';
M = length(v);
N = M - 1;
% fast Chebyshev transforms
b = real(fft([v fliplr(v(2:N))]));
c = real(ifft(1i*[0:M-2 0 2-M:-1].*b));
n2b = (0:M-2).^2.*b(1:N);
% derivative construction
w = zeros(M,1);
w(1) = sum(n2b)/N + 0.5*N*b(M);
w(2:N) = -csc(pi/N*(1:M-2)).*c(2:N);
w(M) = sum((-1).^(1:N).*n2b)/N + 0.5*(-1)^M*N*b(M);

end