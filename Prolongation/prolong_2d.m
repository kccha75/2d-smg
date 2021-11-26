% 2D Restriction operator assuming fine grid is 2x
% coarse grid points
%
% Inputs:
% vc - coarse grid
% x_prolong - prolongation function in x
% y_prolong - prolongation function in y
%
% Outputs:
% vf - coarse grid prolongation


function vf=prolong_2d(vc,x_prolong,y_prolong)

% Prolong in x
vf=x_prolong(vc);

% Prolong in y
vf=y_prolong(vf');

vf=vf';

end