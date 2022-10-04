function w = fchd2(v)
%FCHD2	Fast Chebyshev derivative
% FCHD2(V) computes the first derivative of the data in V located along
% the N+1 Chebyshev points cos(pi*(0:N)/N).
%
% 
% Example 1:
% Use FCHD2 to find the second derivative of the function f(x) = tan(x) 
% over [-1,1], and compare with f''(x) = 2 tan(x) sec(x)^2.
% 
% x = cos(pi*(0:10)/10); % create sparse Chebyshev-spaced grid of 11 points
% xx = linspace(-1,1); % create dense, linearly spaced grid
% plot(xx,2*tan(xx).*sec(xx).^2,x,fchd2(tan(x)));
%
% 
% Example 2:
% To show the spectral convergence property of the Chebyshev derivative,
% compute the error between the Chebyshev second derivative and the exact
% second derivative of f(x) = tan(x) for several N.
% 
% N = 1:30;
% err = zeros(1,length(N));
% 
% for n = N
%     x = cos(pi*(0:n)/n)'; % establish grid
%     err(n) = max(2*tan(x).*sec(x).^2 - fchd2(tan(x))); % compute error
% end
% 
% loglog(N,err); %display
v = v(:)';
M = length(v);
N = M - 1;
% fast Chebyshev transforms
b = real(fft([v fliplr(v(2:N))]));
c = real(ifft(1i*[0:M-2 0 2-M:-1].*b));
d = real(ifft(-[0:M-1 2-M:-1].^2.*b));
% precomputations
theta = pi/N*(1:M-2);
ii2 = (2:N-1).^2;
n2b = ii2.*(ii2-1).*b(3:N)/(N*3);
bM = N*(N^2-1)*b(M)/6;
% derivative construction
w = zeros(M,1);
w(1) = sum(n2b) + bM;
w(2:N) = csc(theta).^2.*(d(2:M-1) - cot(theta).*c(2:M-1));
w(M) = sum((-1).^(2:N-1).*n2b) + (-1)^N*bM;
end