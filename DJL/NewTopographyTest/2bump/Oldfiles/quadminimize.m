% u be a vector ...

function [v,u,i]=quadminimize(V,U,I,DJL,pde,domain,option)

for j=1:50

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

    I(4)=int2d(v,domain,DJL);
    U(4)=u;

    % locate maximum
    [~,kk]=max(abs(I(4)-I));
    % Remove furtheset point ... that isn't the new point
   
    U(kk)=[];
    I(kk)=[];

    % sort ...
    [U,index]=sort([U(1) U(2) U(3)]);
    I=[I(index(1)) I(index(2)) I(index(3))];
    



    
%         % check which end is larger difference
%         if diffU(1)>diffU(end)
%             U(1)=[];
%             I(1)=[];
%         elseif diffU(1)<=diffU(end)
%             U(end)=[];
%             I(end)=[];
%         end
    
%     end

end

end