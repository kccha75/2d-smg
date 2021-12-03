% Alternative code to set boundary conditions on 2D matrix in matrix solve
% 
% Sets boundary conditions after kron (not before)
%
% Code may prove useful if ever extended to mixed derivative terms ...



% Dirichlet conditions, can be easily extended to Neumann

% Top
A(1:Nx:(Ny-1)*Nx+1,:)=sparse(Ny,Nx*Ny); % Zero specific rows to BCs
index=sub2ind(size(A),1:Nx:(Ny-1)*Nx+1,1:Nx:(Ny-1)*Nx+1); % Get index
A(index)=1; % Set main diag element of rows to 1

% Bottom
A(Nx:Nx:Nx*Ny,:)=sparse(Ny,Nx*Ny);
index=sub2ind(size(A),Nx:Nx:Nx*Ny,Nx:Nx:Nx*Ny);
A(index)=1;