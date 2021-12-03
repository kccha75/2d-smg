% Solves au_xx+bu_yy+c*u=f with 2nd-order FD scheme (cheb-fourier mesh)
%
% Uses matrix inversion (backslash)
% Assumes periodic boundary conditions
%
% Inputs:
% v - not used except for size purposes
% pde.a
% pde.b
% pde.f
% domain.dx
% domain.N
% 
% Ouputs:
% v - solution
%

function v=matrixsolve_2d(~,pde,domain,~) % IN PROGRESS

N=domain.N;
k=domain.k;
a=pde.a;
b=pde.b;
c=pde.c;

% Check if coefficients constant
if length(pde.a)==1
    a=pde.a*ones(domain.N);
end
if length(pde.b)==1
    b=pde.b*ones(domain.N);
end
if length(pde.c)==1
    c=pde.c*ones(domain.N);
end

for i=1:length(domain.discretisation)
    
    switch domain.discretisation(i)
        
        case 1 % Fourier
            
            % 1D u_xx
            Dxx{i}=real(ifft(-k{i}.^2.*fft(eye(N(i),N(i)))));
            
            
        case 2 % Cheb
            
            % 1D u_xx
            Dxx{i}=ifct(chebdiff(fct(eye(N(i),N(i))),2));
    
    end
    
    
end

% 2D matrices
D2xx=kron(speye(N(2)),Dxx{1});
D2yy=kron(Dxx{2},speye(N(1)));

% -------------------------------------------------------------------------
% Spectral matrix
A=a(:).*D2xx+b(:).*D2yy+c(:).*speye(Nx*Ny);

% -------------------------------------------------------------------------
% Boundary conditions
% -------------------------------------------------------------------------
% domain.BC will need to have vectors ... and flags of BC type ...

for i=1:length(domain.BCflag)
    
    switch domain.BCflag(i)
        
        
        case 1 % Dirichlet
            
            
        case 2 % Neumann
            
            
        case 3 % Fourier
        
        
        
    end
    
    
end

% Left

% A(1:Nx)=

% Right

% A((Ny-1)*Nx+1:Nx*Ny)=

% Top
A(1:Nx:(Ny-1)*Nx+1,:)=sparse(Ny,Nx*Ny); % Zero specific rows to BCs
index=sub2ind(size(A),1:Nx:(Ny-1)*Nx+1,1:Nx:(Ny-1)*Nx+1); % Get index
A(index)=1; % Set main diag element of rows to 1

% Bottom
A(Nx:Nx:Nx*Ny,:)=sparse(Ny,Nx*Ny);
index=sub2ind(size(A),Nx:Nx:Nx*Ny,Nx:Nx:Nx*Ny);
A(index)=1;

% -------------------------------------------------------------------------
% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.c(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

v=A\pde.f(:);
v=reshape(v,Nx,Ny);

end