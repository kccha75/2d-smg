% Function performs natural parameter continuation
%
% Input:
% v - initial solution at initial parameter
% u - initial parameter
% s - initial step size
% ds - step size
% ds_max - max step size
% ds_min - min step size
% DJL - DJL structure
% domain
% option
%
% Output:
% v - solution vector at each parameter value 
% u - parameter vector

function [v,u]=naturalparametercontinuation(v,u,DJL,domain,cont_option)

steps=cont_option.steps;
ds=cont_option.ds;
ds_min=cont_option.ds_min;
ds_max=cont_option.ds_max;
N_opt=cont_option.N_opt;

% Initialise
u(1)=u;
v(:,:,1)=v;
j=1;
dv=0; % Initial gradient guess ...

while j<steps
    
    % Predictor
    v(:,:,j+1)=v(:,:,j)+dv*ds;

    u(j+1)=u(j)+ds;

    % Update variable
    DJL.u=u(j+1);

    % Update pde / jacobian
    [pde,domain]=DJLpdeinitialise(DJL,domain);

    % Newton iterations
    [v(:,:,j+1),i,flag]=NewtonSolve(v(:,:,j+1),DJL,pde,domain,cont_option);
    
    % Converged and not 0 solution
    if flag==1 && max(max(abs(v(:,:,j+1))))>=1e-8

        fprintf('Converged after %d Newton Iterations step = %d\n',i,j)

        % check overturning
        dv=2*ifct(chebdiff(fct(v(:,:,j+1)'),1));
        if max(dv(:))>1
            fprintf('Overturning detected!\n')
            return
        end

        % Update for next Newton iteration
        dv=(v(:,:,j+1)-v(:,:,j))/ds; % simple dv estimate
        j=j+1;

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

    elseif flag==0 || max(max(abs(v(:,:,j+1))))<1e-8 % not converged or 0 solution
        
        if flag==0
            fprintf('Did not converge to required tolerance  after %d Newton Iterations at step %d\n',i,j)
        end
        if max(max(abs(v(:,:,j+1))))<1e-8
            fprintf('Did not converge to properer solution!!! Converged to 0 after %d Newton Iterations at step %d\n',i,j)
        end
        
        % Halve step size
        ds=ds/2;
        fprintf('Halving step size to %f\n',ds)
        
        % Break loop if minimum step size exceeded
        if ds<=ds_min
            fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
            % Clear last non-converged step
            v(:,:,end)=[];
            u(end)=[];
            return % end function
        end

    end

end