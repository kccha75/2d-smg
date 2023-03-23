
u(1)=u1;
u(2)=DJL.u;

[mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v1,[],[]);
y(1)=maxi-mini;

maxi=v2(maxindex_x,index_z);

mini=v2(minindex_x,index_z);
y(2)=maxi-mini;
i=3;
while abs(y(i-1))>1e-7

    u(i)=u(i-1)-y(i-1)*(u(i-1)-u(i-2))/(y(i-1)-y(i-2));

    DJL.u=u(i);
    v1=v2;
    v2=NewtonSolve(v1,DJL,pde,domain,option);
    [mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v2,minindex_x,maxindex_x);
    y(i)=maxi-mini;
    i=i+1;

end