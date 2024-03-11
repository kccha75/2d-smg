% Function performs Yang's preconditioner for fifth-order KdV equation
% 
% Preconditioner in the form H=a*u_xxxx+b*u_xx+c*u
%
% Chebyshev or Fourier pseudospectral
%
% Uses matrix inversion (backslash)
%
% Inputs:
% v - not used except for size purposes
% pde.a
% pde.b
% pde.f
% domain.dim - can only do 1D or 2D
% domain.discretisation - cheb or fourier
% domain.BC - BC for cheb
% domain.N
% domain.dx
% 
% Ouputs:
% v - solution


function v=yang_kdv_pre(RHS,pde,domain,~)

    k=domain.k{1};

    % Solve!
    if isempty(RHS)
        v=ifft(fft2(pde.f)./(pde.a.*k.^4-pde.b.*k.^2+pde.c));
    else
        v=ifft(fft2(RHS)./(pde.a.*k.^4-pde.b.*k.^2+pde.c));
    end
    
end