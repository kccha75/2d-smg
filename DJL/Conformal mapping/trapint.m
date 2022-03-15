% Function calculates area approximation using trapezoidal rule
%
% area = trapint(f,x)
%
% Inputs:
% f - vector of function f(x) (column vector)
% x - vector of x (column vector)
%
% Outputs:
% area (scalar)

function area=trapint(f,x)

% Formula for trapezoidal rule with summation through matrix multiplcation
area = 0.5*(x(2:end)-x(1:end-1))'*(f(2:end)+f(1:end-1));

end