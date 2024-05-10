% Self contained code for KP (testing purposes)
%

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

% PDE Parameters
a=1;
b=2;
mu=-1.2;
c=-mu;
d=1;

% RHS
f= 0*X;

% Initial guess
v0=-0.43*sech(0.3*sqrt((X).^2+(Y).^2)).*cos(X);


% -------------------------------------------------------------------------
% NEWTON HERE
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=c+6*v0;
v=v0;

% Error guess (keep at 0)
e0=zeros(Nx,Ny);

[KX,KY]=ndgrid(kx,ky);

tol=1e-10;

tic
for i=1:20

    % Initial RHS of linear equation
    F=f-(Lukp(v,a,b,c,d)+3*real(ifft(-kx.^2.*fft(v.^2))));
    
    r=rms(rms(F));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
%     e=gmres(@(u) reshape(Lukp(reshape(u,Nx,Ny),a,b,cnew,d),Nx*Ny,1),F(:),[],tol,1000,@(u) reshape(KPpre(reshape(u,Nx,Ny),a,b,c,d),Nx*Ny,1));
    e=bicgstab(@(u) reshape(Lukp(reshape(u,Nx,Ny),a,b,cnew,d),Nx*Ny,1),F(:),tol,1000,@(u) reshape(KPpre(reshape(u,Nx,Ny),a,b,c,d),Nx*Ny,1));
    e=reshape(e,Nx,Ny);
    
    % Update correction
    v=v+real(e);

    v=fft2(v);
    v(1,1)=0;
    v=real(ifft2(v));

    cnew=c+6*v;

end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc

surf(v);

% -------------------------------------------------------------------------
% Subroutine Lu
% -------------------------------------------------------------------------
function u=Lukp(v,a,b,c,d)
global Nx Ny Lx Ly kx ky KX KY x y X Y mu

u=real(ifft((-a.*kx.^6+b.*kx.^4).*fft(v))-ifft(kx.^2.*fft(c.*v))-transpose(ifft(d.*ky.^2.*fft(transpose(v)))));

end

% -------------------------------------------------------------------------
% Subroutine preconditioner
% -------------------------------------------------------------------------
function u=KPpre(f,a,b,c,d)
global Nx Ny Lx Ly kx ky KX KY x y X Y mu

    C=0.0001;

    u=real(ifft2(fft2(f)./(C+KX.^2.*(a.*KX.^4-b.*KX.^2+c)+1*d.*KY.^2)));

end