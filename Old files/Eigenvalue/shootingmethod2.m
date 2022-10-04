clear;clc;close all

tspan=[0,1];

c=linspace(0,1,10000);

dy=1;

yend=zeros(length(c),1);

for i=1:length(c) % c loop

	y0=[0;1];

	ode=@(z,y) [y(2) ; -z/c(i)^2*y(1)];

	[z,y]=ode45(ode,tspan,y0);
    
	yend(i)=y(end,1);

    
end

plot(c,yend,c,0*c);