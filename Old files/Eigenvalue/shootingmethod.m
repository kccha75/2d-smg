clear;clc;close all

tspan=[0,1];

c=linspace(-1,1,1000);

dy=linspace(-100,100,100);

yend=zeros(length(c),length(dy));

for i=1:length(c) % c loop
    
    for j=1:length(dy) % dy loop
    
        y0=[0;dy(j)];

        ode=@(z,y) [y(2) ; -z/c(i)^2*y(1)];

        [z,y]=ode45(ode,tspan,y0);
    
        yend(i,j)=y(end,1);
        
    end
    
end
surf(c,dy,yend');
xlabel('c');ylabel('dy')
hold on
surf(c,dy,zeros(length(c),length(dy))');