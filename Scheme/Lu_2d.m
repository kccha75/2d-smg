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
function f=Lu_2d(v,pde,domain)

k=domain.k;

a=pde.a;
b=pde.b;
c=pde.c;

u=cell(2);
Lu=u;
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
f=a.*Lu{1}+b.*Lu{2}'+c.*v;

% -------------------------------------------------------------------------
% Apply BCs
% -------------------------------------------------------------------------
% Index cells for each boundary
index=cell(length(domain.BCflag),1);

Nx=domain.N(1);
Ny=domain.N(2);

index{1}=1:Nx:Nx*(Ny-1)+1; % Top boundary of matrix (x(1))
index{2}=Nx:Nx:Nx*Ny; % Bottom boundary of matrix (x(end))
index{3}=1:Nx; % Left boundary of matrix (y(1))
index{4}=Nx*(Ny-1)+1:Nx*Ny; % Right boundary of matrix (y(end))

for i=1:length(domain.BCflag)
    
    switch domain.BCflag(i)
        
        case 1 % Dirichlet
            
            f(index{i})=v(index{i});
            
        case 2 % Neumann
            
            if rem(i,2)~=0 % Is not even?
                f(index{i})=sum(fct(v).*k{i}.^2); % for x=1
            else
                f(index{i})=sum((-1).^k{i}.*fct(v).*k{i}.^2); % for x=-1
            end
            
            
        case 3 % Fourier

            % do nothing :)
    end
    
end

end