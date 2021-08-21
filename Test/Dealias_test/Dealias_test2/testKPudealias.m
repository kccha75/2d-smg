% Test speed ... gives same outputs, similar speed also
n=100;

tic
for i=1:n
z=fourier_KPu_2d_dealias(v,pde,domain);
end
toc
tic
for i=1:n
zz=fourier_KPu_2d_dealias2(v,pde,domain);
end
toc
tic
for i=1:n
zz=fourier_KPu_2d(v,pde,domain);
end
toc