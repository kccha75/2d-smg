% Solves au_xx+bu_yy+c*u=f with 2nd order FD matrix inversion
% OR
% Solves au_xx+bu_x+c*u=f with 2nd order FD matrix inversion
%
% Chebyshev or Fourier pseudospectral
%
% Uses matrix inversion (backslash)
%
% Inputs:
% v - not used except for size purposes
% pde.a
% pde.b
% pde.f
% domain.dim - can only do 1D or 2D
% domain.discretisation - cheb or fourier
% domain.BC - BC for cheb
% domain.N
% domain.dx
% 
% Ouputs:
% v - solution

function v=FDmatrixsolve(~,pde,domain,~)

N=domain.N;
Nx=domain.N(1);
Ny=domain.N(2);
dx=domain.dx;
a=pde.a;
b=pde.b;
c=pde.c;

Dxx=cell(domain.dim,1);
Dx=cell(domain.dim,1);

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
            % Preset diagonal array
            B=zeros(N(i),3);

            B(1:N(i),2)=-2/dx{i}^2;
            B(1:N(i)-1,1)=1/(dx{i}^2);
            B(2:N(i),3)=1/(dx{i}^2);

            % Diagonal position array
            d=[-1 0 1];

            % Generate sparse matrix
            Dxx{i}=spdiags(B,d,N(i),N(i));

            % First row (cyclic)
            Dxx{i}(1,N(i))=1/(dx{i}^2);
            
            % Last row (cyclic)
            Dxx{i}(N(i),1)=1/(dx{i}^2);

            % 1D u_x
            % Preset diagonal array
            B=zeros(N(i),3);

            B(1:N(i),2)=0;
            B(1:N(i)-1,1)=-1/2/(dx{i});
            B(2:N(i),3)=1/2/(dx{i});

            % Diagonal position array
            d=[-1 0 1];

            % Generate sparse matrix
            Dx{i}=spdiags(B,d,N(i),N(i));

            % First row (cyclic)
            Dx{i}(1,N(i))=-1/2/(dx{i});
            
            % Last row (cyclic)
            Dx{i}(N(i),1)=1/2/(dx{i});

        case 2 % Cheb
            
            % 1D u_xx
            % Step size
            dx1=dx{i}(1:end-1);dx2=dx{i}(2:end);

            % Preset diagonal array
            B=zeros(N(i),3);

            % Main diag u_i
            B(2:N(i)-1,2)=-2./(dx1.*dx2);

            % Off diag u_i+1
            B(3:N(i),3)=2./(dx2.*(dx1+dx2));
            
            % Off diag u_i-1
            B(1:N(i)-2,1)=2./(dx1.*(dx1+dx2));
            
            % Diagonal position array
            d=[-1 0 1];

            % Generate sparse matrix
            Dxx{i}=spdiags(B,d,N(i),N(i));

            % 1D u_x
            % Preset diagonal array
            B=zeros(N(i),3);

            % Main diag u_i
            B(2:N(i)-1,2)=(1./dx1-1./dx2);

            % Off diag u_i+1
            B(3:N(i),3)=dx1./(dx2.*(dx1+dx2));
            
            % Off diag u_i-1
            B(1:N(i)-2,1)=-dx2./(dx1.*(dx1+dx2));
            
            % Diagonal position array
            d=[-1 0 1];

            % Generate sparse matrix
            Dx{i}=spdiags(B,d,N(i),N(i));
            
    end
            
end

% -------------------------------------------------------------------------
            % BC matrix
% -------------------------------------------------------------------------
BC_mat=cell(1,domain.dim);
a_mat=cell(3,1);

a_mat{1}=a;
a_mat{2}=b;
a_mat{3}=c;

index=cell(2,domain.dim); % Left and right, for each dimension

index{1,1}=1:Nx:Nx*(Ny-1)+1; % Top boundary of matrix (x(1))
index{2,1}=Nx:Nx:Nx*Ny; % Bottom boundary of matrix (x(end))

index{1,2}=1:Nx; % Left boundary of matrix (y(1))
index{2,2}=Nx*(Ny-1)+1:Nx*Ny; % Right boundary of matrix (y(end))

for i=1:domain.dim
    
    BC_mat{i}=sparse(zeros(N(i)));
    
    if domain.discretisation(i)~=1
        
        % Left neumann x(1) BC
        BC_mat{i}(1,1)=-1/dx{i}(1)-1/(dx{i}(1)+dx{i}(2));
        BC_mat{i}(1,2)=1/dx{i}(1)+1/dx{i}(2);
        BC_mat{i}(1,3)=-dx{i}(1)/(dx{i}(2)*(dx{i}(1)+dx{i}(2)));
        BC_mat{i}(1,:)=BC_mat{i}(1,:).*domain.BC{2,i}; % BC coefficient
        % (forward difference scheme i,i+1,i+2)

        % Dirichlet 
        BC_mat{i}(1,1)=BC_mat{i}(1,1)+domain.BC{1,i};

        % Right neumann x(end) BC
        BC_mat{i}(end,end)=1/dx{i}(end)+1/(dx{i}(end)+dx{i}(end-1));
        BC_mat{i}(end,end-1)=-1/dx{i}(end)-1/dx{i}(end-1);
        BC_mat{i}(end,end-2)=dx{i}(end)/(dx{i}(end-1)*(dx{i}(end)+dx{i}(end-1)));
        BC_mat{i}(end,:)=BC_mat{i}(end,:).*domain.BC{4,i}; % BC coefficient  
        % (backward difference scheme i,i-1,i-2)

        % Dirichlet 
        BC_mat{i}(end,end)=BC_mat{i}(end,end)+domain.BC{3,i};
        
        for j=1:3
            a_mat{j}(index{1,i})=0;
            a_mat{j}(index{2,i})=0;
        end
        
    end
    
end
% -------------------------------------------------------------------------
% 1D matrices
if domain.dim==1

    % Spectral matrix
    A=a_mat{1}(:).*Dxx{1}+a_mat{2}(:).*Dx{1}+a_mat{3}(:).*speye(Nx*Ny)+BC_mat{1};

end
% -------------------------------------------------------------------------
% 2D matrices
if domain.dim==2

    % BC matrix
    BC2_mat=kron(speye(N(2)),BC_mat{1})+kron(BC_mat{2},speye(N(1)));

    D2xx=kron(speye(N(2)),Dxx{1});
    D2yy=kron(Dxx{2},speye(N(1)));

    % Spectral matrix
    A=a_mat{1}(:).*D2xx+a_mat{2}(:).*D2yy+a_mat{3}(:).*speye(Nx*Ny)+BC2_mat;

end
% -------------------------------------------------------------------------
% Solve!

v=A\pde.f(:);
v=reshape(v,Nx,Ny);

end