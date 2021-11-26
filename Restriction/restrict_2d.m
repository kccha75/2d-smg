% 2D Restriction operator assuming fine grid is 2x
% coarse grid points
%
% Inputs:
% vf - fine grid
% x_restrict - restriction function in x
% y_restrict - restriction function in y
%
% Outputs:
% vc - fine grid restriction


function vc=restrict_2d(vf,x_restrict,y_restrict)

% Restrict in x
vc=x_restrict(vf);

% Restrict in y
vc=y_restrict(vc');
vc=vc';


end