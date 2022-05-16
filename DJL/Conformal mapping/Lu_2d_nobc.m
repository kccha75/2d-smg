% Performs multiplication L*u using Spectral collocation methods
% L is the spectral operator for the given PDE in the form
% au_xx+bu_yy+cu=f
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% domain.k - wave number
%
% Ouputs:
% f - L*u
%
function f=Lu_2d_nobc(v,pde,domain)

k=domain.k;

a=pde.a;
b=pde.b;
c=pde.c;

u=cell(2);
Lu=cell(2);
u{1}=v;
u{2}=v';

for i=1:domain.dim
    
    switch domain.discretisation(i)
        
        case 1 % Fourier
            
            % 1D u_xx
            Lu{i}=real(ifft(-k{i}.^2.*fft(u{i})));
            
        case 2 % Cheb
            
            % 1D u_xx
            Lu{i}=ifct(chebdiff(fct(u{i}),2));
    end
    
    
end

% Compute au_xx+bu_yy+cu matrix
f=a.*Lu{1}+b.*Lu{2}'+c.*u{1};

end