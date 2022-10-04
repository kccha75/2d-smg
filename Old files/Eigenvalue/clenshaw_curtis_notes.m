% Input f - vector at chebyshev points


function I = clenshaw_curtis(f)
n=length(f)-1;
g = real(fft(f([1:n+1 n:-1:2])/(2*n)));          % FFT
a = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)];  % Chebyshev coefficients
w = 0*a';
w(1:2:end) = 2./(1-(0:2:n).^2);  % weight vector
    I = w*a;                                   % the integral