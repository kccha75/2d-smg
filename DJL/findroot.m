% Function finds index of location where v(index)~=x0
%
% Inputs: 
% v - array
% x0 - desired value
%
% Outputs:
% index - where v(index)~=x0
% such that x0 is between v(index) and v(index+1)

function index=findroot(v,x0)

% Turn into sign change problem
v=v-x0;

% Detect sign change index
index = find((v(1:end-1)>0 & v(2:end)<0) | (v(1:end-1)<0 & v(2:end)>0));

end