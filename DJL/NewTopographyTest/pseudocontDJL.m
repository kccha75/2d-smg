% Function performs pseudo arclength continuation on DJL equation
% au_xx+bu_zz+N^2(z-u)u/c^2=0
%
% Updates parameter c in continuation 
%
% Input:
% v - initial solution at initial parameter
% dv - initial solution gradient
% lambda - initial parameter
% dlambda - initial parameter gradient

% domain
% option - solver options
% cont_option - continuation options
%
% Output:
% V - solution vector at each parameter value 
% U - parameter vector

function [V,U]=pseudocontDJL(v,dv,lambda,dlambda,DJL,domain,option,cont_option)

ds=cont_option.ds; % Step size initial
ds_min=cont_option.ds_min; % Minimum step size
ds_max=cont_option.ds_max; % Maximum step size
N_opt=cont_option.N_opt; % Optimum Newton iterations per step
Newtonmaxit=cont_option.Newtonmaxit; % Maximum Newton iterations per step
steps=cont_option.steps; % maximum steps allowed

Newtontol=option.Newtontol;
tailtol=option.tailtol;

% Set max Newton iterations to continuation option
option.Newtonmaxit=Newtonmaxit;

% Initialise
V(:,:,1)=v;
U(1)=lambda;

j=1;

% Newton parameters
flag=0;
e0=zeros(domain.N); % Error guess (keep at 0)

while j<steps
    
    % Predictor
    V(:,:,j+1)=V(:,:,j)+dv*ds;

    U(j+1)=U(j)+dlambda*ds;

    % Update variable
    DJL.u=U(j+1);

    % Update pde / jacobian
    [DJL,pde,domain]=DJLpdeinitialise_topography(DJL,domain);

% -------------------------------------------------------------------------
% Newton here
% -------------------------------------------------------------------------
    for i=1:Newtonmaxit

        % Calculate Jacobian for linear equation
        J=option.jacobian(V(:,:,j+1),DJL,pde,domain);
        
        % Check nonlinear residual
        r=rms(J.f(:));
        fprintf('Residual Newton = %d\n',r)

        % Solve linear equation
        RHS1=J.f;
        RHS2=-2*DJL.N2((domain.X{2}+1)/2-V(:,:,j+1)).*V(:,:,j+1)/U(j+1)^3; RHS2(:,end)=0;
        
        J.f=RHS1;
        z1=mg(e0,J,domain,option);
%         z1=bicgstab(e0,J,domain,option);
        J.f=RHS2;
        z2=mg(e0,J,domain,option);
%         z2=bicgstab(e0,J,domain,option);

        if r<=Newtontol

            flag=1;
            numNewtonit=i-1;
            break

        end
        
        % Update correction
        delta_lambda=(ds-sum(dot(dv,(V(:,:,j+1)-V(:,:,j))))-dlambda*(U(j+1)-U(j))-sum(dot(dv,z1))) ...
            /(dlambda-sum(dot(dv,z2)));
        delta_v=z1-delta_lambda*z2;
    
        V(:,:,j+1)=V(:,:,j+1)+delta_v;
        U(j+1)=U(j+1)+delta_lambda;

        % Update variable
        DJL.u=U(j+1);

        % Update pde / jacobian
        [DJL,pde,domain]=DJLpdeinitialise_topography(DJL,domain);

    end

    % If converged final loop ...
    % Calculate Jacobian for linear equation
    J=option.jacobian(V(:,:,j+1),DJL,pde,domain);

    if rms(J.f(:))>Newtontol
    
        fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
        flag=0;

    elseif flag==0 && i==Newtonmaxit

        fprintf('Converged after %d Newton Iterations \n',i)
        flag=1;
        numNewtonit=i;
    
    end 

% -------------------------------------------------------------------------

    % Converged newton and not 0 solution and decays to 0
    if flag==1 && max(max(abs(V(:,:,j+1))))>tailtol && max(max(abs(V(1,:,j+1))))<tailtol
        
        fprintf('Converged after %d Newton Iterations step = %d\n',numNewtonit,j)

        % check overturning
        diffv=2*real(ifct(chebdiff(real(fct(transpose(V(:,:,j+1)))),1)));
        if max(diffv(:))>1
            fprintf('Overturning detected!\n')
            return
        end

        % Update for next Newton iteration
        dlambda=1/(dlambda-sum(dot(dv,z2)));
        dv=-dlambda*z2;
    
        % Normalise
        mag=sqrt(sum(dot(dv,dv))+dlambda^2);
        dlambda=dlambda/mag;
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

    elseif flag==0 || max(max(abs(V(:,:,j+1))))<=tailtol || max(max(abs(V(1,:,j+1))))>=tailtol
        
        if flag==0
            fprintf('Did not converge to required tolerance after %d Newton Iterations at step %d\n',i,j)
        end

        if abs(max(V(1,:,j+1)))>=1e-8
            fprintf('Non decaying solution at step %d\n',j)
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

end