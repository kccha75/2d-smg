% Function uses Preconditioned BiConjugate Gradient Stability to solve system Au=f
%
% Iterates until tolerance condition OR maximum iterations reached
%
% Inputs:
% v - best estimate
% pde - structure consisting of pde coefficients
% scheme.k - wave number
% option.operator - operator used to find Lu
% option.numit - number of iterations
% option.preconditioner - preconditioner used
% option.prenumit - number of preconditioner iterations
%
% Ouputs:
% v - solution after tol / maxit reached
% r - residual for solution
%
% Optional display messages can be commented out or left in
%
% NOTE: only left preconditioning at the moment as input, can take left and
% right preconditioning in algorithm

function [v,r]=bicgstab(v,pde,domain,option)
f=pde.f;
maxit=100;

% Initial residual
r=pde.f-option.operator(v,pde,domain);
p=r;
d=r;

if rms(r(:))<option.tol
    fprintf('Initial guess already below tolerance!\n');
    return
end


for i=1:maxit
    
    % Check if preconditioner available (note single use of preconditioner)
    % should be K1*K2y=p generally where K=K1*K2
    if ~isempty(option.preconditioner) && option.prenumit~=0
    
        pde.f=p;
        option.numit=option.prenumit;
        y=option.preconditioner(zeros(size(v)),pde,domain,option);
    
    else
    
        y=p;
    
    end

    delta=dot(r(:),d(:));
    
    q=option.operator(y,pde,domain);
    alpha=delta/dot(d(:),q(:));
    s=r-alpha*q;
    
    % Tolerance check
    if rms(s(:))<option.tol
        v=v+alpha*p;
        fprintf('Bi Conjugate Gradient STAB Converged after %d iterations!\n',i);
        break
    end
    
    
% Check if preconditioner available (note single use of preconditioner)
% should be K1*K2z=s generally where K=K1*K2
if ~isempty(option.preconditioner) && option.prenumit~=0
    
    pde.f=s;
    option.numit=option.prenumit;
    z=option.preconditioner(zeros(size(v)),pde,domain,option);
    
    t=option.operator(z,pde,domain);% A*s
    
    pde.f=t;
    K1t=option.preconditioner(t,pde,domain,option); % inv(K1)*t
    pde.f=s;
    K1s=option.preconditioner(s,pde,domain,option); % inv(K1)*s
    
    w=dot(K1t(:),K1s(:))/dot(K1t(:),K1t(:));
    
else
    
    z=s;
    
    t=option.operator(z,pde,domain);% A*s
    w=dot(t(:),s(:))/dot(t(:),t(:));
    
end 
    
    v=v+alpha*y+w*z; % estimate of new solution
    
    % Recalculate every 100 iteration to remove roundoff errors
    if mod(i,20)==0
        
        r=f-option.operator(v,pde,domain);
        
    else
        
        r=s-w*t;%disp(rms(r(:)));
        
    end
    
    % Tolerance check
    if rms(r(:))<option.tol
         
        fprintf('Bi Conjugate Gradient STAB Converged after %d iterations!\n',i);
        break
        
    end
     
    deltanew=dot(r(:),d(:));
    beta=(alpha/w)*deltanew/delta;
    p=r+beta*(p-w*q);
    
%     % Check r.d =/= 0
%     if abs(dot(r(:),d(:)))<option.tol
% %         disp(abs(dot(r(:),d(:))))
%         d=r;
%         p=r;
%     end
    
end

if i==maxit
    
    fprintf('Bi Conjugate Gradient STAB did not converge after %d iterations!\n',i);
    fprintf('Residual %d\n',rms(rms(r)));
    
end

end