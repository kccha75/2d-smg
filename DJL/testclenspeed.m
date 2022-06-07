N=2^5;

k=(0:N)';
x=cos(pi*k/(N));
f=exp(x);

tic
[x,w]=clencurt(N);
I1=w*f;
toc

tic
I2=clenshaw_curtis(f);
toc