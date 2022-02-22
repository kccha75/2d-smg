% Performs multiplication L*u using Spectral collocation methods
% L is the spectral operator for the given PDE in the form
% au_xx+bu=f
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% domain.k - wave number
%
% Ouputs:
% f - L*u
%
function f=Lu_1d(v,pde,domain)

k=domain.k{1};

a=pde.a;
b=pde.b;

switch domain.discretisation
        
    case 1 % Fourier
            
        % 1D u_xx
        Lu=real(ifft(-k.^2.*fft(v)));
            
        
    case 2 % Cheb
            
        % 1D u_xx
        Lu=ifct(chebdiff(fct(v),2));
    
end


% Compute au_xx+bu matrix
f=a.*Lu+b.*v;

% -------------------------------------------------------------------------
% Apply BCs
% -------------------------------------------------------------------------
if domain.discretisation~=1
	f(1)=domain.BC{1,1}.*v(1)+ ...
        domain.BC{2,1}.*sum(fct(v).*k.^2); % x(1)
	f(end)=domain.BC{3,1}.*v(end)+ ... 
        domain.BC{4,1}.*sum((-1).^k.*fct(v).*k.^2); % x(end)
end


end