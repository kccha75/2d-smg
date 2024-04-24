clear
global Nx Ny Lx Ly kx ky KX KY x y X Y mu
Nx=2^10;
Ny=2^7;

Lx=120*pi;
Ly=60*pi;

kx=2*pi/Lx*[0:Nx/2-1 -Nx/2 -Nx/2+1:-1]';
x=Lx*(-Nx/2:Nx/2-1)'/Nx;

ky=2*pi/Ly*[0:Ny/2-1 -Ny/2 -Ny/2+1:-1]';
y=Ly*(-Ny/2:Ny/2-1)'/Ny;

[X,Y] = ndgrid(x,y);

mu=-1.2;

% Initial guess
v0=-0.43*sech(0.3*sqrt((X).^2+(Y).^2)).*cos(X);


% -------------------------------------------------------------------------
% NEWTON HERE
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=mu+6*v0;
v=v0;

% Error guess (keep at 0)
e0=zeros(Nx,Ny);

[KX,KY]=ndgrid(kx,ky);

tol=1e-10;

tic
for i=1:20

    % Initial RHS of linear equation
    F=-(Lukp(v)+3*real(ifft2(-KX.^2.*fft2(v.^2))));
    
    r=rms(rms(F));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end

    % Solve linear equation
%     e=gmres(@(u) reshape(Lukp(reshape(u,Nx,Ny),a,b,cnew,d),Nx*Ny,1),F(:),[],tol,1000,@(u) reshape(KPpre(reshape(u,Nx,Ny),a,b,c,d),Nx*Ny,1));
    e=bicgstab(@(u) reshape(Lukp(reshape(u,Nx,Ny)),Nx*Ny,1)+6*reshape(real(ifft2(-KX.^2.*fft2(v))),Nx*Ny,1).*u,F(:),tol,1000,@(u) reshape(KPpre(reshape(u,Nx,Ny)),Nx*Ny,1));
    e=reshape(e,Nx,Ny);
    
    % Update correction
    v=v+real(e);

    v=fft2(v);
    v(1,1)=0;
    v=real(ifft2(v));

end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc

surf(v);

% -------------------------------------------------------------------------
% Subroutine Lu
% -------------------------------------------------------------------------
function u=Lukp(v)
global Nx Ny Lx Ly kx ky KX KY x y X Y mu

u=real(ifft2(fft2(v).*(-KX.^2.*(KX.^4-2*KX.^2-mu)-KY.^2)));

end

% -------------------------------------------------------------------------
% Subroutine preconditioner
% -------------------------------------------------------------------------
function u=KPpre(f,a,b,c,d)
global Nx Ny Lx Ly kx ky KX KY x y X Y mu

    C=0.0001;

    u=real(ifft2(fft2(f)./(C+KX.^2.*(KX.^4-2*KX.^2-mu)+0*KY.^2)));

end