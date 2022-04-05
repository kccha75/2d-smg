% Function performs natural parameter continuation
%
% Input:
% v - initial solution at initial parameter
% u - initial parameter
% s - initial step size
% ds - step size
% dsmax
% dsmin
% DJL
% domain
% option
%
% Output:
% v - solution vector at each parameter value 
% u - parameter vector

% function [v,lambda]=naturalparametercontinuation(v,u,du,domain)
steps=cont_option.steps;
ds=cont_option.ds;
ds_min=cont_option.ds_min;
ds_max=cont_option.ds_max;
N_opt=cont_option.N_opt;

N=domain.N;
x=domain.x;
dx=domain.dx;

u(1)=pde.u;

j=1;
v(:,:,1)=v;

while j<steps
    
    % Predictor
    if j==1
        v(:,:,j+1)=v(:,:,j);
    else
        v(:,:,j+1)=v(:,:,j)+(v(:,:,j)-v(:,:,j-1));
    end

    u(j+1)=u(j)+ds;

    % Update variable
    pde.u=u(j+1);
    DJL.u=u(j+1);

    % Update pde / jacobian
    [pde,domain]=DJL_pde_initialise(DJL,domain);

    % Newton iterations
    [v(:,:,j+1),i,flag]=NewtonSolve(v(:,:,j+1),pde,domain,cont_option);
    
    % Converged and not 0 solution
    if flag==1 && max(abs(v(:)))>=1e-10

        fprintf('Converged after %d Newton Iterations step = %d\n',i,j)

        % check overturning
        dv=2*ifct(chebdiff(fct(v(:,:,j+1)'),1));
        if max(dv(:))>1
            fprintf('Overturning detected!')
            return
        end

        % Optimal step length control for next step
            xi=N_opt/i;
            
            if xi<0.5
                ds=ds*0.5;
            elseif xi>=0.5 && xi<=2
                ds=ds*xi;
            elseif xi>2
                ds=ds*2;
            end
            
            % Max step size check
            if ds>ds_max
                ds=ds_max;
            elseif ds<ds_min
                fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
                return
            end
            fprintf('New step size to %f\n',ds)

            % Update variable for next Newton iteration
            j=j+1;
            v(:,:,j+1)=v(:,:,j);
            u(j+1)=u(j)+ds;
        
            % Update variable
            pde.u=u(j+1);
            DJL.u=u(j+1);

            % Update pde / jacobian
            [pde,domain]=DJL_pde_initialise(DJL,domain);

    elseif flag==0 || max(abs(v(:)))<1e-10 % or 0 solution
        
        if flag==0
            fprintf('Did not converge to required tolerance  after %d Newton Iterations at step %d\n',i,j)
        end
        if max(abs(v(:)))<1e-10
            fprintf('Did not converge to properer solution!!! Converged to 0 after %d Newton Iterations at step %d\n',i,j)
            pause
        end
        % Halve step size
        ds=ds/2;
        fprintf('Halving step size to %f\n',ds)
        
        % Break loop if minimum step size exceeded
        if ds<=ds_min
            fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
            return % end function
        end

    end


end

