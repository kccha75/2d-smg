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

function v=matrixsolve_2d(~,pde,domain,~)

N=domain.N;
Nx=domain.N(1);
Ny=domain.N(2);
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
            
% -------------------------------------------------------------------------            
            % APPLY BCs
% -------------------------------------------------------------------------  

            if i==1
                
                a(1,:)=1;a(end,:)=1;
                b(1,:)=0;b(end,:)=0;
                c(1,:)=0;c(end,:)=0;
                
            elseif i==2
                
                a(1,:)=0;a(end,:)=0;
                b(1,:)=1;b(end,:)=1;
                c(1,:)=0;c(end,:)=0;
                
            end

            for j=2*i:-1:2*i-1
                
                switch domain.BCflag(j)
                    
                    case 1 % Dirichlet
                        
                        if mod(j,2)==1 % x(1) boundary
                        	Dxx{i}(1,:)=0;
                            Dxx{i}(1,1)=1;
                        else % x(end) boundary
                            Dxx{i}(end,:)=0;
                            Dxx{i}(end,end)=1;
                        end
                        
                    case 2 % Neumann
                        
                        if mod(j,2)==1 % x(1) boundary
                            Dxx{i}(1,:)=sum(fct(eye(N(i),N(i))).*k{i}.^2);
                        else            % x(end) boundary
                            Dxx{i}(end,:)=sum((-1).^k{i}.*fct(eye(N(i),N(i))).*k{i}.^2);
                        end
                    
                end
                
            end
    
    end
    
    
end

% -------------------------------------------------------------------------
% 2D matrices
D2xx=kron(speye(N(2)),Dxx{1});
D2yy=kron(Dxx{2},speye(N(1)));

% Spectral matrix
A=a(:).*D2xx+b(:).*D2yy+c(:).*speye(Nx*Ny);
% -------------------------------------------------------------------------
% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.c(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

v=A\pde.f(:);
v=reshape(v,Nx,Ny);

end