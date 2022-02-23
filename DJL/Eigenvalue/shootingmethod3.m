clear;clc;close all

tspan=[0,1];

lambda=linspace(0,1000,10000);

dy=1;

yend=zeros(length(lambda),1);

for i=1:length(lambda) % c loop

	y0=[0;1];

	ode=@(z,y) [y(2) ; -z*lambda(i)*y(1)];

	[z,y]=ode45(ode,tspan,y0);
    
	yend(i)=y(end,1);

    
end

plot(lambda,yend,lambda,0*lambda);