% Function solves nonlinear pde using Newton's method
%
% Inputs:
% 
% v0 - initial guess
% nonlinear - nonlinear functions for the jacobian and nonlinear residual
% pde - pde coefficients
% domain.N 
% option.jacobian - jacobian function
% option.numit=1 - iterations
% option.tol - linear residual tolerance
% option.Newtontol - Newton set tolerance
% option.Newtonmaxit - Newton max iterations
% option.Newtonsolver - Linear solver for Newton iterations
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

    if r>1e10 || isnan(r)
        flag=0;
        numNewtonit=i-1;
        return
    end

    if r<=option.Newtontol

        fprintf('Converged after %d Newton Iterations \n',i-1)
        flag=1;
        numNewtonit=i-1;
        return

    end
    
    % Solve linear equation
%     option.tol=1e-3*r;

    switch option.Newtonsolver
        case 'mg'

            [e,r]=mg(e0,J,domain,option);

        case 'gmres'
            
            func=@(v) reshape(Lu(reshape(v,domain.N),J,domain),prod(domain.N),1);
            pre=@(v) reshape(FDmatrixsolve(v,J,domain,[]),prod(domain.N),1);
            [e,flag,r,iter]=gmres(func,J.f(:),[],option.tol,option.numit,pre);
            e=reshape(e,domain.N);
            
        case 'bicgstab'

            func=@(v) reshape(Lu(reshape(v,domain.N),J,domain),prod(domain.N),1);
            pre=@(v) reshape(FDmatrixsolve(v(:),J,domain,[]),prod(domain.N),1);
            [e,flag,r,iter]=bicgstab(func,J.f(:),option.tol,option.numit,pre);
            e=reshape(e,domain.N);

        otherwise

            fprintf('Invalid Newton solver!\n')
            return
    end
   
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
