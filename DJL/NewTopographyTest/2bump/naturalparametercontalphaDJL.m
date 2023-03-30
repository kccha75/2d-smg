% Function performs natural parameter continuation
%
% Input:
% v1 - initial solution at initial parameter
% u1 - initial parameter
% v2 - initial solution at initial parameter 2
% u2 - initial parameter 2
% cont_option.ds - step size
% cont_option.ds_max - max step size
% cont_option.ds_min - min step size
% DJL - DJL structure
% domain
% option - solver option
%
% Output:
% V - solution vector at each parameter value 
% U - parameter vector

function [V,U,W]=naturalparametercontalphaDJL(v,delta,alpha,DJL,domain,option,cont_option)

ds=cont_option.ds;
ds_min=cont_option.ds_min;
ds_max=cont_option.ds_max;
N_opt=cont_option.N_opt;
Newtonmaxit=cont_option.Newtonmaxit;
steps=cont_option.steps;

Newtontol=option.Newtontol;
tailtol=option.tailtol;

% Set max Newton iterations to continuation option
option.Newtonmaxit=Newtonmaxit;

% Initialise
V(:,:,1)=v;
W(1)=delta;
U(1)=alpha;

j=1;
dv=0; % Initial gradient guess ...
dw=0;

while j<steps
    
    % Predictor
    V(:,:,j+1)=V(:,:,j)+dv*ds;
    W(j+1)=W(j)+dw*ds;
    U(j+1)=U(j)+ds;

    % Update variable
    DJL.alpha=U(j+1);
    DJL.u=W(j+1);
    DJL.v=V(:,:,j+1);

    % Update mapping
    [DJL,domain]=conformalmapping(DJL,domain,option);

    % Update pde / jacobian
    [DJL,pde,domain]=DJLpdeinitialise_topography(DJL,domain);

    % Newton 1 here ...
    [v1,i,newtonflag1]=NewtonSolve(DJL.v,DJL,pde,domain,option);
    
    % Newton 2 small delta adjustment
    DJL.u=DJL.u-1e-4;
    [v2,~,newtonflag2]=NewtonSolve(v1,DJL,pde,domain,option);

    % Update secant 
    [V(:,:,j+1),W(j+1),y,secantflag]=DJLtabletopsecant(v1,v2,W(j+1),W(j+1)-1e-4,DJL,pde,domain,option);

    % Check secant convergence
    if secantflag==0
            fprintf('Secant did not converge\n')
    end

    % Converged and not 0 solution
    if newtonflag1==1 && newtonflag2==1 && secantflag==1 &&...
            max(max(abs(v1)))>tailtol && max(max(abs(v2)))>tailtol

        fprintf('Converged after %d Newton Iterations step = %d\n',i,j)

        % check overturning
        diffv=2*real(ifct(chebdiff(real(fct(transpose(V(:,:,j+1)))),1)));
        if max(diffv(:))>1
            fprintf('Overturning detected!\n')
            return
        end

        % Update for next Newton iteration
        dv=(V(:,:,j+1)-V(:,:,j))/ds; % simple dv estimate
        dw=(W(j+1)-W(j))/ds; % simple dw estimate
        j=j+1;

        % Optimal step length control for next step (minus first..)
        if j>3
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
        end

    elseif newtonflag1==0 || newtonflag2==0 || secantflag==0 || max(max(abs(V(:,:,j+1))))<=tailtol || max(max(abs(V(1,:,j+1))))>=tailtol
        
        if newtonflag1==0 || newtonflag2==0

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
            W(end)=[];
            return % end function
        end

    end

end