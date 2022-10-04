% Function generates preconditioning matrix of operator A using second
% order finite differencing

% second order FD on -u_xx+au_x+bu

% Inputs:
% L - domain size
% N - grid size
% a - a(x) array
% b - b(x) array

% Outputs:
% H - preconditioning matrix (sparse)

function H=Hfe(L,N,a,b)

% Step size
dx=L/N;

B(1:N,2)=4;
B(1:N-1,1)=1;
B(2:N,3)=1;

% Diagonal position array
d=[-1 0 1];

% Generate sparse matrix
M=spdiags(B,d,N,N);

% First row (cyclic)
M(1,N)=1;

% Last row (cyclic)
M(N,1)=1;

M=M*1/6*dx;


K=Hfd(L,N,a,b)*dx;

H=inv(M)*K;

end