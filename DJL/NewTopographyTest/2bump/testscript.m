[U,index]=sort([u u1 u2]);
V(:,:,1)=v;
V(:,:,2)=v2;
V(:,:,3)=v1;
I=[int2d(v.^2,domain,DJL) int2d(v2.^2,domain,DJL) int2d(v1.^2,domain,DJL)];

[v,u,i]=quadminimize(V,U,I,DJL,pde,domain,option);