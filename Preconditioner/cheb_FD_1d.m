% Solves -u_xx+a*u_x+b*u=f with 2nd-order FD scheme (NON UNIFORM MESH)
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

function v=cheb_FD_1d(~,pde,domain,~)

N=domain.N(1);
dx=domain.dx{1};
a=-pde.a;
b=-pde.b;

% Check if coefficients constant
if length(pde.a)==1
    a=a*ones(domain.N);
end
if length(pde.b)==1
    b=b*ones(domain.N);
end

% Step size
dx1=dx(1:end-1);dx2=dx(2:end);

% Preset diagonal array
B=zeros(N,3);

% Main diag u_i
B(2:N-1,2)=(1./dx1+1./dx2).*2./(dx1+dx2)-dx1.*a(2:N-1)./(dx2.*(dx2+dx1)) ...
    +dx2.*a(2:N-1)./(dx1.*(dx1+dx2))+b(2:N-1);

% Off diag u_i+1
B(3:N,3)=-2./(dx2.*(dx1+dx2))+a(2:N-1).*dx1./(dx2.*(dx2+dx1));

% Off diag u_i-1
B(1:N-2,1)=-2./(dx1.*(dx1+dx2))-a(2:N-1).*dx2./(dx1.*(dx2+dx1));

% Diagonal position array
d=[-1 0 1];

% Generate sparse matrix
A=spdiags(B,d,N,N);

% Dirichlet BCs
A(1,1)=1;A(end,end)=1;

% Solve
v=-A\pde.f;

end