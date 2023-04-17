% u be a vector ...

function [v,u,i]=quadminimize(V,U,I,DJL,pde,domain,option)

for j=1:5

    % Quadratic interpolation ...
    
    p=polyfit(U,I,2);
    
    % Turning point of quadratic
    u=-p(2)/(2*p(1));
    
    % Find solution at this point ...
    % initial guess the closest one ...
    [~,k]=min(abs(u-U));
    v=V(:,:,k);
    
    % Check if minimum value ...
    if min(u-U<0)==1 && min(abs(u-U))>max(abs(diff(U)))
        u=U(k)-max(abs(diff(U)));
    elseif min(u-U>0)==1 && min(abs(u-U))>max(abs(diff(U)))
        u=U(k)+max(abs(diff(U)));
    end

    % Newton ...
    DJL.u=u;
    [v,numit,newtonflag]=NewtonSolve(v,DJL,pde,domain,option);

    % Check Newton converge
    if newtonflag==0
        fprintf('AHHHHHH\n')
        return
    end
    i=int2d(v,domain,DJL);

    % Remove furtheset point ... that isn't the new point
    [U,index]=sort([U(1) U(2) U(3) u]);
    I=[I(1) I(1) I(3) i];
    diffU=diff(U);
    
    % u is smallest
    if index(1)==4
        % remove largest
        U(end)=[];
        I(end)=[];
    
    % u is largest
    elseif index(4)==4 
        % remove smallest
        U(1)=[];
        I(1)=[];
    
    else
    
        % check which end is larger difference
        if diffU(1)>diffU(end)
            U(1)=[];
            I(1)=[];
        elseif diffU(1)<=diffU(end)
            U(end)=[];
            I(end)=[];
        end
    
    end

end

end