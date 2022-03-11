% -------------------------------------------------------------------------
% NEWTON HERE
% -------------------------------------------------------------------------

v=v0;

% Error guess (keep at 0)
e0=zeros(N);

tic
for i=1:10
    
    % Calculate Jacobian for linear equation
    jacobian=jacobian_DJL(v,pde,domain);
   
    % Check nonlinear residual
    r=rms(rms(jacobian.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
%     option.tol=1e-3*r;
    [e,r]=bicgstab(e0,jacobian,domain,option);
%     [e,r]=mg(e0,jacobian,domain,option);

    % Update correction
    v=v+e;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc

