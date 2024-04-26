% Function uses Preconditioned Bi-Conjugate Gradient STAB to solve system 
% Au=f
%
% Iterates until tolerance condition OR maximum iterations reached
%
% Inputs:
% v - best estimate
% pde - structure consisting of pde coefficients
% scheme.k - wave number in x
% option.operator - operator used to find Lu
% option.numit - number of iterations
% option.preconditioner - preconditioner used
% option.prenumit - number of preconditioner iterations
%
% Ouputs:
% v - solution after tol / maxit reached
% r - residual for solution
%
% Optional display messages can be commented out or left in
%
% Uses MATLAB inbuilt function

function [v,r]=bicgstabmg(v,pde,domain,option)
c=pde.c;
A=@(u) reshape(option.operator(reshape(u,domain.N),pde,domain),prod(domain.N),1);
pde.c=1.2;
M=@(u) reshape(option.preconditioner(reshape(u,domain.N),pde,domain),prod(domain.N),1);

[v,flag,r]=bicgstab(A,pde.f(:),option.tol,100,M);

v=reshape(v,domain.N);
pde.c=c;
end