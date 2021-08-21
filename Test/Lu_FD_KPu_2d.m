% Solves (au_xx+bu)_xx+du_yy with 2nd-order FD scheme (midpoints)
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

function Lu=Lu_FD_KPu_2d(v,pde,domain)

Nx=domain.N(1);
Ny=domain.N(2);
dx=domain.dx(1);
dy=domain.dx(2);

% 4th derivative 2nd order FD matrix in x
B=zeros(Nx,5);
B(:,1)=1;
B(:,2)=-4;
B(:,3)=6;
B(:,4)=-4;
B(:,5)=1;

D4x=spdiags(B,[-2,-1,0,1,2],Nx,Nx);

% Periodic BC
% First row cyclic
D4x(1,Nx)=-4;
D4x(1,Nx-1)=1;
% Second row
D4x(2,Nx)=1;
% Last row
D4x(Nx,1)=-4;
D4x(Nx,2)=1;
% Second last row
D4x(Nx-1,1)=1;

D4x=D4x/dx^4;

% 2th derivative 2nd order FD matrix in y
B=zeros(Ny,3);
B(:,1)=1;
B(:,2)=-2;
B(:,3)=1;

D2y=spdiags(B,[-1,0,1],Ny,Ny);

% Periodic BC
% First row cyclic
D2y(1,Ny)=1;
% Last row
D2y(Ny,1)=1;

D2y=D2y/dy^2;

% 2th derivative 2nd order FD matrix in x
B=zeros(Nx,3);
B(:,1)=1;
B(:,2)=-2;
B(:,3)=1;

D2x=spdiags(B,[-1,0,1],Nx,Nx);

% Periodic BC
% First row cyclic
D2x(1,Nx)=1;
% Last row
D2x(Nx,1)=1;

D2x=D2x/dx^2;

%(au_xx+bu)_xx+du_yy
A=pde.a(:).*kron(speye(Ny),D4x)+pde.b(:).*kron(speye(Ny),D2x)+pde.d(:).*kron(D2y,speye(Nx));

Lu=A*v(:)+pde.c(:).*kron(speye(Ny),D2x)*(v(:).^2);
Lu=reshape(Lu,Nx,Ny);
end