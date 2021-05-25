y2=cos(pi*(0:Nx-1)/(Nx-1))';

y2c= cos(pi*(0:Nx/2-1)/(Nx/2-1))';

x2=Lx*(-Nx/2:Nx/2-1)/Nx;

x2c=x2(1:2:end);

[X2,Y2]=ndgrid(x2,y2);

u=sin(X2).*(Y2+1).*(Y2-1);

[X2C,Y2C]=ndgrid(x2c,y2c);

uc=sin(X2C).*(Y2C+1).*(Y2C-1);

uu=fourier_cheb_restrict(u);

surf(x2c,y2c,uu);figure;surf(x2c,y2c,uc)
