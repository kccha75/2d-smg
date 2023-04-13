

function objectparameter=objectivefunction(v,u,DJL,pde,domain,option)

if u>.36 || u < .34
    u=NaN;
    objectparameter=NaN;
    return
end

DJL.u=u;
[v,i,flag]=NewtonSolve(v,DJL,pde,domain,option);

if flag ==0

    fprintf('Newton did not converge ...\n')
    u=NaN;

end

objectparameter=int2d(v.^2,domain,DJL);

end