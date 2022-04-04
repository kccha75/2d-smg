% Function solves nonlinear pde using Newton's method
%
% Inputs:
% 
% v0 - initial guess
% pde - pde coefficients
% domain.N 
% option.jacobian - jacobian function
% option.Newtontol - Newton set tolerance
% option.maxit - Newton max iterations
%
% Outputs:
%
% v - solution after Newton iterations
% i - iterations required
% flag - 1 if converged to tolerance
%      - 0 if after max iteration did not reach tolerance  

function [v,numNewtonit,flag]=NewtonSolve(v0,pde,domain,option)

v=v0;

% Error guess (keep at 0)
e0=zeros(domain.N);

tic
for i=1:option.maxit
    
    % Calculate Jacobian for linear equation
    jacobian=option.jacobian(v,pde,domain);
   
    % Check nonlinear residual
    r=rms(jacobian.f(:));
    fprintf('Residual Newton = %d\n',r)

    if r<=option.Newtontol

        fprintf('Converged after %d Newton Iterations \n',i-1)
        flag=1;
        break

    end
    
    % Solve linear equation
%     option.tol=1e-3*r;
    [e,r]=bicgstab(e0,jacobian,domain,option);
%     [e,r]=mg(e0,jacobian,domain,option);

    % Update correction
    v=v+e;
    
end

if i==option.maxit
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    flag=0;
    
end 

numNewtonit=i;

end
