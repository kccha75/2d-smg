% Function performs natural parameter continuation on epsilon in N2
% (transition from linear to sech2)
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

function [V,U,W]=naturalparametercontalphaDJLepsilon(v,delta,epsilon,DJL,domain,option,cont_option)
global psi0 d
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
U(1)=epsilon;

j=1;
dv1=0; % Initial gradient guess ...
dw1=0;
ds1=ds;
skip=1;
while j<steps
    
    if j<=2
        % Predictor
        V(:,:,j+1)=V(:,:,j)+dv1*ds;
    else
        V(:,:,j+1)=V(:,:,j)+dv2*ds;
    end
    
%     W(j+1)=W(j)+dw*ds;
    U(j+1)=U(j)+ds;

    % Update variable
    DJL.mu=U(j+1);
%     DJL.u=W(j+1);
    DJL.v=V(:,:,j+1);

    DJL.Lx=2*DJL.KAI/DJL.mu^2;

    % Update mapping
    [DJL,domain]=conformalmapping(DJL,domain,option);

    epsilon=U(j+1);
    DJL.N2=@(psi) (epsilon*(-psi+1)+(1-epsilon)*sech((psi-psi0)/d).^2)/(epsilon/2-d*(epsilon-1)*(tanh((1-psi0)/d)+tanh(psi0/d)));
    DJL.N2d=@(psi) (epsilon*(-1)+(1-epsilon)*(-2*sech((psi-psi0)/d).^2.*tanh((psi-psi0)/d))/d)/(epsilon/2-d*(epsilon-1)*(tanh((1-psi0)/d)+tanh(psi0/d)));

    % Update pde / jacobian
    [DJL,pde,domain]=DJLpdeinitialise_topography(DJL,domain);

    % Newton 1 here ...
    DJL.u=W(j)+dw1*ds;
    [v1,~,newtonflag1]=NewtonSolve(DJL.v,DJL,pde,domain,option); % important for BCs
    
    if j<=2 % at first step 
        % Newton 2 small delta adjustment
        DJL.u=W(j)+1e-5;
        [v2,newtonit,newtonflag2]=NewtonSolve(v1,DJL,pde,domain,option);

        if newtonflag1==1 && newtonflag2==1
            % Update secant 
            [V(:,:,j+1),W(j+1),y,i,secantflag]=DJLtabletopsecant(v1,v2,W(j)+dw1*ds,W(j)+1e-5,DJL,pde,domain,option);
        else
            secantflag=0;
            i=0;
        end
    
    else

        % Quadratic extrapolation?
        DJL.u=W(j)+dw2*ds;
        [v2,newtonit,newtonflag2]=NewtonSolve(v1,DJL,pde,domain,option);

        % if v1=v2 ...
        if newtonit==0
            DJL.u=W(j)+1e-5;
            [v2,~,newtonflag2]=NewtonSolve(v1,DJL,pde,domain,option);
        end

        if newtonflag1==1 && newtonflag2==1
            % Update secant 
            [V(:,:,j+1),W(j+1),y,i,secantflag]=DJLtabletopsecant(v1,v2,W(j)+dw1*ds,W(j)+dw2*ds,DJL,pde,domain,option);
        else
            secantflag=0;
            i=0;
        end

    end

    % Check secant convergence
    if secantflag==0
            fprintf('Secant did not converge\n')
    end

    % Converged and not 0 solution
    if secantflag==1 &&...
            max(max(abs(v1)))>tailtol && max(max(abs(v2)))>tailtol

        fprintf('Converged after %d Newton Iterations step = %d\n',i,j)

%         check overturning
%         diffv=2*real(ifct(chebdiff(real(fct(transpose(V(:,:,j+1)))),1)));
%         if max(diffv(:))>1
%             fprintf('Overturning detected!\n')
%             return
%         end

        % Update for next Newton iteration
        dv1=(V(:,:,j+1)-V(:,:,j))/ds; % simple dv estimate
        dw1=(W(j+1)-W(j))/ds; % simple dw estimate
        
        if j~=1

            ds2=ds1; % ds2 is the distance between j and j-1
            ds1=ds; % ds1 is the distance between j+1 and j
%             dv2=(3/2*V(:,:,j+1)-2*V(:,:,j)+1/2*V(:,:,j-1))/ds;
%             dw2=(3/2*W(j+1)-2*W(j)+1/2*W(j-1))/ds;

            dv2=(V(:,:,j-1)*ds1^2-V(:,:,j)*(ds1+ds2)^2+V(:,:,j+1)*ds2*(2*ds1+ds2))/(ds1*ds2*(ds1+ds2));
            dw2=(W(j-1)*ds1^2-W(j)*(ds1+ds2)^2+W(j+1)*ds2*(2*ds1+ds2))/(ds1*ds2*(ds1+ds2));

            % Minimum 'step size' for w
            if abs(dw1-dw2)<1e-5
                dw2=dw1+sign(dw1-dw2)*1e-5;
            end
        end
        skip=1;
        j=j+1;
        
        % Optimal step length control for next step (minus first..)
        if j>2
            xi=N_opt/i;
            
            if xi<0.5
                ds=ds*0.5;
            elseif xi>=0.5 && xi<=2
                ds=ds*xi;
            elseif xi>2
                ds=ds*2;
            end
            
            % Max step size check
            if abs(ds)>ds_max
                ds=sign(ds)*ds_max;
            elseif abs(ds)<ds_min

%                 fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
%                 return
            end
            fprintf('New step size to %f\n',ds)
        end
    
    elseif secantflag==0 || max(max(abs(V(:,:,j+1))))<=tailtol || max(max(abs(V(1,:,j+1))))>=tailtol

        if secantflag==0

            fprintf('Secant did not converge\n')

        end

        
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
        if abs(ds)<=ds_min
            if skip==1
                ds=sign(ds)*0.1;
                skip=0;
            else

                fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
                % Clear last non-converged step
                V(:,:,end)=[];
                U(end)=[];
                W(end)=[];
                return % end function
            end
        end

    end

end