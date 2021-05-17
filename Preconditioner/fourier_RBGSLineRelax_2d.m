% Performs RBGS Line relaxation on -(a*u_x)_x-(b*u_y)_y+c*u=f
% with 2nd-order FD approximation
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% pde.f
% domain.dx
% domain.N
% option.numit - number of iterations
%
% Ouputs:
% v - best estimate after relaxation
%
% Notes:
% Assumes constant dx dy (for now)

function v=fourier_RBGSLineRelax_2d(v,pde,domain,option)

for i=1:option.numit
    
    % Line Relax y direction
    v  = fourier_RBGSLineRelaxY(v,pde,domain);
    
    % Line Relax x direction (swap a,b + transpose f)
    temp=pde.a;
    pde.a=pde.b;
    pde.b=temp;
    pde.f=pde.f';
    v  = fourier_RBGSLineRelaxY(v',pde,domain);
    
    % Transpose for solution
    v = v';

end

end

% -------------------------------------------------------------------------
% Subroutine Line Relax Y direction
% -------------------------------------------------------------------------

function v = fourier_RBGSLineRelaxY(v,pde,domain)

Nx=domain.N(1);
Ny=domain.N(2);
dx=domain.dx(1);
dy=domain.dx(2);

index1x=[(2:Nx) 1];
index2x=(1:Nx);
index3x=[Nx (1:Nx-1)];

% Check if coefficients constant and find midpoint of a and b (averaging)
if length(pde.a)==1
    
    a_mid=pde.a*ones(Nx,Ny);
    
else
    
    a_mid=(pde.a(index1x,:)+pde.a(index2x,:))/2;
    
end

index1y=[(2:Ny) 1];
index2y=(1:Ny);

if length(pde.b)==1
    
    b_mid=pde.b*ones(Nx,Ny);
    
else
    
    b_mid=(pde.b(:,index1y)+pde.b(:,index2y))/2;
    
end

g=zeros(Nx,Ny);

for n=1:2
    
    col = n:2:Ny; % i
    col1 = mod(col,Ny)+1; % i+1
    col2 = mod(col-2,Ny)+1; % i-1
    
    % Generate RHS
    g(:,col)=pde.f(:,col)+1/dy^2*(b_mid(:,col).*v(:,col1)+b_mid(:,col2).*v(:,col2));
    
    % Generate LHS tridiagonal matrix + solve
    for i=1:length(col)
        
        B=zeros(Nx,3);
        B(index2x,2)=(a_mid(index2x,col(i))+a_mid(index3x,col(i)))/dx^2+(b_mid(:,col(i))+b_mid(:,col1(i)))/dy^2+pde.c(:,col(i)); % main diag
        B(1:Nx-1,1)=-a_mid(1:Nx-1,col(i))/dx^2; % lower diag
        B(2:Nx,3)=-a_mid(1:Nx-1,col(i))/dx^2; % upper diag
    
        A=spdiags(B,-1:1,Nx,Nx);
        
        % Periodic condition
        A(1,Nx)=-a_mid(Nx,col(i))/dx^2;
        A(Nx,1)=-a_mid(Nx,col(i))/dx^2;
    
        v(:,col(i))=A\g(:,col(i));
    end
    
end

end