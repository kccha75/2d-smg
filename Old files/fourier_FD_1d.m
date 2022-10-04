% Solves -u_xx+a*u_x+b*u=f with 2nd-order FD scheme
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
% Notes:
% Assumes constant dx dy (for now)

function v=fourier_FD_1d(~,pde,domain,~)

N=domain.N(1);
dx=domain.dx(1);
a=pde.a;
b=pde.b;

% Check if coefficients constant
if length(pde.a)==1
    a=a*ones(domain.N);
end
if length(pde.b)==1
    b=b*ones(domain.N);
end

% Preset diagonal array
B=zeros(N,3);

B(1:N,2)=(2+b(1:N)*dx^2)/dx^2;
B(1:N-1,1)=(-2-a(2:N)*dx)/(2*dx^2);
B(2:N,3)=(-2+a(1:N-1)*dx)/(2*dx^2);

% Diagonal position array
d=[-1 0 1];

% Generate sparse matrix
A=spdiags(B,d,N,N);

% First row (cyclic)
A(1,N)=(-2-a(1)*dx)/(2*dx^2);

% Last row (cyclic)
A(N,1)=(-2+a(N)*dx)/(2*dx^2);

% Solve
v=A\pde.f;

end