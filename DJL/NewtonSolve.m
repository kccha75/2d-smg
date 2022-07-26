% Function solves nonlinear pde using Newton's method
%
% Inputs:
% 
% v0 - initial guess
% nonlinear - nonlinear functions for the jacobian and nonlinear residual
% pde - pde coefficients
% domain.N 
% option.jacobian - jacobian function
% option.tol - linear residual tolerance
% option.Newtontol - Newton set tolerance
% option.Newtonmaxit - Newton max iterations
%
% Outputs:
%
% v - solution after Newton iterations
% numNewtonit - iterations required
% flag - 1 if converged to tolerance
%      - 0 if after max iteration did not reach tolerance  

function [v,numNewtonit,flag]=NewtonSolve(v,nonlinear,pde,domain,option)

flag=0;

% Error guess (keep at 0)
e0=zeros(domain.N);

tic
for i=1:option.Newtonmaxit
    
    % Calculate Jacobian for linear equation
    J=option.jacobian(v,nonlinear,pde,domain);
   
    % Check nonlinear residual
    r=rms(J.f(:));
    fprintf('Residual Newton = %d\n',r)

    if r<=option.Newtontol

        fprintf('Converged after %d Newton Iterations \n',i-1)
        flag=1;
        numNewtonit=i-1;
        return

    end
    
    % Solve linear equation
%     option.tol=1e-3*r;
    [e,r]=bicgstab(e0,J,domain,option);
%     [e,r]=mg(e0,J,domain,option);

    % Update correction
    v=v+e;
    
end

% If converged final loop ...
% Calculate Jacobian for linear equation
J=option.jacobian(v,nonlinear,pde,domain);

if rms(J.f(:))>option.Newtontol
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    flag=0;

elseif flag==0

    fprintf('Converged after %d Newton Iterations \n',i)
    flag=1;
    
end 

numNewtonit=i;

end
