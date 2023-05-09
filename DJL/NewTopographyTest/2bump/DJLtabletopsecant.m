% Function locates tabletop solution location using secant method
%
% Assumes max and min are at the same location, this difference can be
% negative. Uses old solution as new guess
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
% i - iterations number
% flag - 1 if converged to tolerance
%      - 0 if after max iteration did not reach tolerance  

function [v,u,y,i,flag]=DJLtabletopsecant(v1,v2,u1,u2,DJL,pde,domain,option)

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
flag=0;
i=1;
j=0;
secant=1;

% incase loop 1 coverge
v=v2;
u=u2;

% loop until difference small, or divide by 0 due to 0 Newton iterations
while i<=option.Newtonmaxit

    % Check convergence ..
    if (abs(y(i+1))<option.Newtontol)

        fprintf('Converged after %d Secant Iterations \n',i-1)
        flag=1;
        break

    elseif y(i+1)-y(i)==0 && abs(y(i+1))<1e-4

        fprintf('0 Newton Iterations detected!\n')
        flag=1;
        break

    end

    if secant==1
        % Secant method
        u=u2-y(i+1)*(u2-u1)/(y(i+1)-y(i));
    end
    % Update delta
    DJL.u=u;

    % Solve for new solution
    [v,numit,newtonflag]=NewtonSolve(v2,DJL,pde,domain,option);

    % Check Newton converge
    if newtonflag==0
        fprintf('Newton did not converge in secant method!\n')

        % smaller step!
        u=(u+u2)/2;
        fprintf('Reducing u to %d\n',u)
        secant=0;
        if abs(u-u2)<abs(u1-u2)
            fprintf('Secant fail!\n')
            flag=0;
            return
        end

    elseif newtonflag==1

        % Find max-min difference
        [mini,maxi,minindex_x,maxindex_x,index_z]=locatemaxmin(v,minindex_x,maxindex_x);
        y(i+2)=maxi-mini;
        fprintf('Secant difference is %d\n',y(i+2))

        % Next loop update
        u1=u2;
        u2=u;
        v2=v;
        i=i+1;
        secant=1;

    end
end

% incase last loop
% Check convergence ..
if (abs(y(end))<option.tol) && i==option.Newtonmaxit

    fprintf('Converged after %d Secant Iterations \n',option.Newtonmaxit)
    flag=1;

elseif y(end)-y(end-1)==0 && i==option.Newtonmaxit && abs(y(end))<1e-4

    fprintf('0 Newton Iterations detected!\n')
    flag=1;

end

end