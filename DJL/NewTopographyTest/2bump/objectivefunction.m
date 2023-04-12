

function objectparameter=objectivefunction(v,u,DJL,pde,domain,option)

if u>.37 || u < .32
    u=NaN;
    objectparameter=NaN;
    return
end

DJL.u=u;
[v,i,flag]=NewtonSolve(v,DJL,pde,domain,option);

if flag ==0

    fprintf('Newton did not converge ...\n')
    return

end

objectparameter=int2d(v,domain,DJL);

end