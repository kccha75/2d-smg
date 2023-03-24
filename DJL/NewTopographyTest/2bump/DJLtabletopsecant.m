% Function locates and 
%
% Inputs:
% v1 - DJL solution 1
% v2 - DJL solution 2
% u1 - delta at DJL solution 1
% u2 - delta at DJL solution 2
% DJL - in Newton
% pde - in Newton
% domain - in Newton
% option - in Newton
%
% Outputs:
% v - DJL tabletop solution (hopefully)
% u - delta value of DJL tabletop solution
% y - difference vector at each step

function [v,u,y]=DJLtabletopsecant(v1,v2,u1,u2,DJL,pde,domain,option)

% vector of delta
u(1)=u1;
u(2)=u2;

% find local minimum and local maximum
[mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v1,[],[]);

% Difference vector (max - min) at location (note can be negative)
y(1)=maxi-mini;

% use index_x to find max and min for v2
maxi=v2(maxindex_x,index_z);
mini=v2(minindex_x,index_z);

% Difference vector (max - min) at location (note can be negative)
y(2)=maxi-mini;

% initialise loop
i=3;
v=v1;
% loop until difference small, or divide by 0 due to 0 Newton iterations
while (abs(y(i-1))>option.tol) && (y(i-1)-y(i-2)~=0)

    % Secant method
    u(i)=u(i-1)-y(i-1)*(u(i-1)-u(i-2))/(y(i-1)-y(i-2));
    
    % Update delta
    DJL.u=u(i);

    % Solve for new solution
    v=NewtonSolve(v,DJL,pde,domain,option);

    % Find max-min difference
    [mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v,minindex_x,maxindex_x);
    y(i)=maxi-mini;
    fprintf('Difference is %d\n',y(i))
    i=i+1;

end

u=u(end);

end