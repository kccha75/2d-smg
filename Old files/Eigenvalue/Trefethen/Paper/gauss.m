function I = gauss(f,n)% (n+1)-pt Gauss quadrature of f
beta = .5./sqrt(1-(2*(1:n)).^(-2)); % 3-term recurrence coeffs

T = diag(beta,1) + diag(beta,-1);   % Jacobi matrix
[V,D] = eig(T); % eigenvalue decomposition
x = diag(D); [x,i] = sort(x);   % nodes (= Legendre points)
w = 2*V(1,i).^2;    % weights
I = w*feval(f,x);   % the integral
