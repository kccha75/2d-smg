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
% V - solution vector at each parameter value 
% U - parameter vector

function [V,U]=pseudocont(v,dv,u,du,fKdV,domain,cont_option)

topography=fKdV.topography;
X=domain.X;

Newtontol=cont_option.Newtontol;
Newtonmaxit=cont_option.Newtonmaxit;
steps=cont_option.steps;
ds=cont_option.ds;
ds_min=cont_option.ds_min;
ds_max=cont_option.ds_max;
N_opt=cont_option.N_opt;


% Initialise
U(1)=u;
V(:,1)=v;

j=1;

% Newton parameters
flag=0;
e0=zeros(domain.N); % Error guess (keep at 0)

while j<steps
    
    % Predictor
    V(:,j+1)=V(:,j)+dv*ds;

    U(j+1)=U(j)+du*ds;

    % Update variable
    fKdV.u=U(j+1);

    % Update pde / jacobian
    [pde,domain]=fKdVpdeinitialise(fKdV,domain);

% -------------------------------------------------------------------------
% Newton here
% -------------------------------------------------------------------------
    for i=1:Newtonmaxit

        % Calculate Jacobian for linear equation
        J=cont_option.jacobian(v,fKdV,pde,domain);
        
        % Check nonlinear residual
        r=rms(J.f(:));
        fprintf('Residual Newton = %d\n',r)

        if r<=Newtontol

            fprintf('Converged after %d Newton Iterations \n',i-1)
            flag=1;
            numNewtonit=i-1;
            break

        end

        % Solve linear equation
        RHS1=-J.f;
        RHS2=-topography(X);

        J.f=RHS1;
        z1=mg(e0,J,domain,cont_option);
        J.f=RHS2;
        z2=mg(e0,J,domain,cont_option);
        
        % Update correction
        delta_u=(ds-dot(dv,(V(:,j+1)-V(:,j)))-du*(U(j+1)-U(j))-dot(dv,z1)) ...
            /(du-dot(dv,z2));
        delta_v=z1-delta_u*z2; 
    
        v=v+delta_v;
        u=u+delta_u;

    end

    % If converged final loop ...
    % Calculate Jacobian for linear equation
    J=cont_option.jacobian(v,fKdV,pde,domain);

    if rms(J.f(:))>option.Newtontol
    
        fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
        flag=0;

    elseif flag==0

        fprintf('Converged after %d Newton Iterations \n',i)
        flag=1;
    
    end 

        numNewtonit=i;

% -------------------------------------------------------------------------

    % Converged newton
    if flag==1
        
        fprintf('Converged after %d Newton Iterations step = %d\n',i,j)

        % Update for next Newton iteration
        du=1/(du-dot(dv,z2));
        dv=-du*z2;
    
        % Normalise
        mag=sqrt(dot(dv,dv)+du^2);
        du=du/mag;
        dv=dv/mag;

        j=j+1;

        % Optimal step length control for next step
            xi=N_opt/numNewtonit;
            
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

    elseif flag==0
        
        if flag==0
            fprintf('Did not converge to required tolerance  after %d Newton Iterations at step %d\n',i,j)
        end
        
        % Halve step size
        ds=ds/2;
        fprintf('Halving step size to %f\n',ds)
        
        % Break loop if minimum step size exceeded
        if ds<=ds_min
            fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
            % Clear last non-converged step
            V(:,:,end)=[];
            U(end)=[];
            return % end function
        end

    end

end