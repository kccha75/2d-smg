% Performs multiplication L*u
% Chebyshev or Fourier pseudospectral methods
%
% L is the spectral operator for the given PDE in the form
%
% au_xx+bu_yy+c*u=f
% OR
% au_xx+bu_x+c*u=f 
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% domain.dim
% domain.discretisation
% domain.BC
% domain.k - wave number
%
% Ouputs:
% f - L*u

function f=Lu_nobc(v,pde,domain)

k=domain.k;

a=pde.a;
b=pde.b;
c=pde.c;

u=cell(domain.dim,1);
Dx=cell(domain.dim,1);
Dxx=cell(domain.dim,1);

u{1}=v;
if domain.dim==2;u{2}=v';end

for i=1:domain.dim
    
    switch domain.discretisation(i)
        
        case 1 % Fourier

            % 1D u_xx
            Dxx{i}=real(ifft(-k{i}.^2.*fft(u{i})));

            % 1D u_x
            Dx{i}=real(ifft(1i*k{i}.*fft(u{i})));

        case 2 % Cheb

            % 1D u_xx
            Dxx{i}=real(ifct(chebdiff(real(fct(u{i})),2)));
            % 1D u_x
            Dx{i}=real(ifct(chebdiff(real(fct(u{i})),1)));

    end
    
end

if domain.dim==1
    % Compute au_xx+bu_x+cu matrix
    f=a.*Dxx{1}+b.*Dx{1}+c.*u{1};
end

if domain.dim==2
    % Compute au_xx+bu_yy+cu matrix
    f=a.*Dxx{1}+b.*Dxx{2}'+c.*u{1};
end


end