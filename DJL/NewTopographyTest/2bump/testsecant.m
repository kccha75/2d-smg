u(1)=u1;
u(2)=DJL.u;

y(1)=locatemaxmin(v1);
y(2)=locatemaxmin(v2);

i=3;
while y(i-1)>1e-5

    u(i)=u(i-1)-y(i-1)*(u(i-1)-u(i-2))/(y(i-1)-y(i-2));

    DJL.u=u(i);
    v1=v2;
    v2=NewtonSolve(v1,DJL,pde,domain,option);
    y(i)=locatemaxmin(v2);

    i=i+1;

end