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

function f=Lu(v,pde,domain)

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
            Dxx{i}=ifct(chebdiff(fct(u{i}),2));
            % 1D u_x
            Dx{i}=ifct(chebdiff(fct(u{i}),1));

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

% -------------------------------------------------------------------------
% Apply BCs
% -------------------------------------------------------------------------
index=cell(2,domain.dim); % Left and right, for each dimension

Nx=domain.N(1);
Ny=domain.N(2);

index{1,1}=1:Nx:Nx*(Ny-1)+1; % Top boundary of matrix (x(1))
index{2,1}=Nx:Nx:Nx*Ny; % Bottom boundary of matrix (x(end))

index{1,2}=1:Nx; % Left boundary of matrix (y(1))
index{2,2}=Nx*(Ny-1)+1:Nx*Ny; % Right boundary of matrix (y(end))

for i=1:domain.dim
    
    if domain.discretisation(i)~=1
        f(index{1,i})=domain.BC{1,i}.*v(index{1,i})+ ...
            domain.BC{2,i}.*sum(fct(u{i}).*k{i}.^2); % x(1)
        f(index{2,i})=domain.BC{3,i}.*v(index{2,i})+ ... 
            domain.BC{4,i}.*sum((-1).^(k{i}+1).*fct(u{i}).*k{i}.^2); % x(end)
    end
    
end

end