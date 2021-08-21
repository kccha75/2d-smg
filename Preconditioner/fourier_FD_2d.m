% Solves -(a*u_x)_x-(b*u_y)_y+c*u=f with 2nd-order FD scheme (midpoints)
%
% Uses matrix inversion (backslash)
% Assumes periodic boundary conditions
%
% Inputs:
% v - not used except for size purposes
% pde.a
% pde.b
% pde.c
% pde.f
% domain.dx
% domain.N
% 
% Ouputs:
% v - solution
%
% Notes:
% Assumes constant dx dy (for now)

function v=fourier_FD_2d(~,pde,domain,~)

Nx=domain.N(1);
Ny=domain.N(2);
dx=domain.dx(1);
dy=domain.dx(2);

% Check if coefficients constant and find midpoint of a and b (at i+1)
if length(pde.a)==1
    
    a_mid=pde.a*ones(Nx,Ny);
    
else

    a_mid=(pde.a+circshift(pde.a,[0,-1]))/2;
    
end

if length(pde.b)==1
    
    b_mid=pde.b*ones(Nx,Ny);
    
else
    
    b_mid=(pde.b+circshift(pde.b,[-1,0]))/2;
    
end

% 1D Fourier Dx diff matrix at i+1/2
B=zeros(Nx,2);
B(:,1)=-1/dx;
B(:,2)=1/dx;
D1xf=spdiags(B,[0,1],Nx,Nx);
D1xf(Nx,1)=B(Nx,2); % periodic condition

% 2D Fourier Dx diff matrix at i+1/2
D2xf=kron(speye(Ny),D1xf);

% 1D Fourier Dy diff matrix at i+1/2
B(:,1)=-1/dy;
B(:,2)=1/dy;
D1yf=spdiags(B,[0,1],Ny,Ny);
D1yf(Ny,1)=B(Ny,2); % periodic condition

% 2D Fourier Dy diff matrix at i+1/2
D2yf=kron(D1yf,speye(Nx));

% 1D Fourier Dx diff matrix at i-1/2
B(:,1)=-1/dx;
B(:,2)=1/dx;
D1xb=spdiags(B,[-1,0],Nx,Nx);
D1xb(1,Nx)=B(1,1); % periodic condition

% 2D Fourier Dx diff matrix at i-1/2
D2xb=kron(speye(Ny),D1xb);

% 1D Fourier Dy diff matrix at i-1/2
B(:,1)=-1/dy;
B(:,2)=1/dy;
D1yb=spdiags(B,[-1,0],Ny,Ny);
D1yb(1,Ny)=B(1,1); % periodic condition

% 2D Fourier Dy diff matrix at i-1/2
D2yb=kron(D1yb,speye(Nx));

% Solve + reshape
A=-D2xb*(a_mid(:).*D2xf)-D2yb*(b_mid(:).*D2yf)+pde.c(:).*speye(Nx*Ny);

% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.c(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

v=A\pde.f(:);
v=reshape(v,Nx,Ny);

end