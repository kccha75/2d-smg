% Solves -u_xx+a*u_x+b*u=f using Fourier matrix inversion
%
% Inputs:
% pde.a
% pde.b
% pde.f
% domain.N
% domain.k
% 
% Ouputs:
% v - solution

function v=fourier_matrixsolve_1d(~,pde,domain,~)

N=domain.N(1);
kx=domain.k{1};

% Fourier 1st derivative matrix
D1x=ifft(1i.*kx.*fft(eye(N,N)));

% Fourier 2nd derivative matrix
D2x=ifft(-kx.^2.*fft(eye(N,N)));

% Operator diff matrix
A=-D2x+pde.a.*D1x+pde.b.*eye(N,N);

% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.b(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

% Solve
v=real(A\pde.f);

end