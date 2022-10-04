Nx=32;

x = cos(pi*(0:Nx-1)/(Nx-1))';

xc= cos(pi*(0:Nx/2-1)/(Nx/2-1))';

xf = cos(pi*(0:2*Nx-1)/(2*Nx-1))';

y=(x-1).*(x+1).^3;

yc=cheb_restrict(y);

yf=cheb_prolong(y);

plot(x,y,'x-',xc,yc,'o-',xf,yf,'*-')


