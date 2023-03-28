% Function locates tabletop solution location using secant method
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
% flag - 1 if converged to tolerance
%      - 0 if after max iteration did not reach tolerance  

function [v,u,y,flag]=DJLtabletopsecant(v1,v2,u1,u2,DJL,pde,domain,option)

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
v=v1;
flag=0;

% loop until difference small, or divide by 0 due to 0 Newton iterations
for i=3:option.Newtonmaxit

    % Check convergence ..
    if (abs(y(i-1))<option.tol)

        fprintf('Converged after %d Secant Iterations \n',i-1)
        flag=1;
        break

    elseif y(i-1)-y(i-2)==0

        fprintf('0 Newton Iterations detected!\n')
        flag=1;
        break

    end

    % Secant method
    u(i)=u(i-1)-y(i-1)*(u(i-1)-u(i-2))/(y(i-1)-y(i-2));
    
    % Update delta
    DJL.u=u(i);

    % Solve for new solution
    [v,numit,newtonflag]=NewtonSolve(v,DJL,pde,domain,option);

    % Check Newton converge
    if newtonflag==0
        fprintf('Newton did not converge in secant method!\n')
        flag=0;
        u=u(end);
        return;
    end

    % Find max-min difference
    [mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v,minindex_x,maxindex_x);
    y(i)=maxi-mini;
    fprintf('Secant difference is %d\n',y(i))

end

% incase last loop
if numit==0
    fprintf('0 Newton Iterations detected!\n')
    flag=1;
end

if (abs(y(end))<option.tol)

    fprintf('Converged after %d Secant Iterations \n',i-1)
    flag=1;
end
u=u(end);

end