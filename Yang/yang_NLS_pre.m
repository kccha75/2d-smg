% Function performs Yang's preconditioner for NLS
% 
% Preconditioner in the form H=c-au_xx-bu_yy
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


function v=yang_NLS_pre(~,pde,domain,~)

    c=3;
    [KX,KY]=ndgrid(domain.k{1},domain.k{2});
    K2=pde.a.*KX.^2+pde.b.*KY.^2;
    v=ifft2(fft2(pde.f)./(c+K2));
    
end